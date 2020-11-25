#include "RapidExternalEvtGen.h"

#include <cstdlib>
#include <fstream>
#include <queue>

#include "TRandom.h"
#include "TSystem.h"

#ifdef RAPID_EVTGEN
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"

#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

bool RapidExternalEvtGen::decay(std::vector<RapidParticle*>& parts) {
#ifdef RAPID_EVTGEN
	EvtParticle* theParent(0);

	if(parts.size() < 1) {
		std::cout << "WARNING in RapidExternalEvtGen::decay : There are no particles to decay." << std::endl;
		return false;
	}
	EvtId theId = EvtPDL::evtIdFromLundKC(parts[0]->id());

	// Parent particle 4-momentum
	TLorentzVector pIn = parts[0]->getP();
	EvtVector4R pInit(pIn.E(),pIn.Px(),pIn.Py(),pIn.Pz());

	theParent = EvtParticleFactory::particleFactory(theId, pInit);
	if (theParent->getSpinStates() == 3) {theParent->setVectorSpinDensity();}

	// Generate the event
	evtGen_->generateDecay(theParent);

	// Store particles to read in the order RapidSim stores them
	std::queue<EvtParticle*> evtParts;
	// Also store the number of children expected for each of these particles so we can remove PHOTOS photons
	std::queue<int> nExpectedChildren;
	evtParts.push(theParent);
	nExpectedChildren.push(parts[0]->nDaughters());

	EvtVector4R x4Evt;
	EvtVector4R p4Evt;
	TLorentzVector p4TLV;

	int iPart=1; // The momentum and origin vertex of the first particle are already set

	while(!evtParts.empty()) {
		EvtParticle* theParticle = evtParts.front();

		// B0 and Bs may mix in EvtGen - RapidSim ignores this step and only records the second state
		while(theParticle->getNDaug()==1) {
			theParticle = theParticle->getDaug(0);
		}

		uint nChildren = nExpectedChildren.front();

		// Loop over the daughter tracks
		for (uint iChild = 0; iChild < nChildren; ++iChild) {
			EvtParticle* child = theParticle->getDaug(iChild);

			if (child != 0) {
				p4Evt = child->getP4Lab();
				x4Evt = child->get4Pos();
				p4TLV.SetPxPyPzE(p4Evt.get(1),p4Evt.get(2),p4Evt.get(3),p4Evt.get(0));
				if(parts.size() < iPart+1u) {
					std::cout << "WARNING in RapidExternalEvtGen::decay : EvtGen has produced too many particles." << std::endl;
					return false;
				}
				parts[iPart]->setP(p4TLV);
				parts[iPart]->getOriginVertex()->setXYZ(x4Evt.get(1),x4Evt.get(2),x4Evt.get(3));
				evtParts.push(child);
				nExpectedChildren.push(parts[iPart]->nDaughters());
				++iPart;
			}
		}

		// Clean up any PHOTOS photons
		for (uint iChild = nChildren; iChild < theParticle->getNDaug(); ++iChild) {
			delete theParticle->getDaug(iChild);
		}
		delete theParticle;
		evtParts.pop();
		nExpectedChildren.pop();
	}

	return true;
#else
	if(!suppressWarning_) {
		std::cout << "WARNING in RapidExternalEvtGen::decay : EvtGen extension not compiled. Will not use EvtGen to decay " << parts[0]->name() << "." << std::endl;
		suppressWarning_=true;
	}

	return false;
#endif
}

bool RapidExternalEvtGen::setup() {
#ifdef RAPID_EVTGEN
	std::cout << "INFO in RapidExternalEvtGen::setup : Setting decay for external EvtGen generator." << std::endl;
	if(!evtGen_) setupGenerator();
	evtGen_->readUDecay(decFileName_.Data());
	return true;
#else
	std::cout << "WARNING in RapidExternalEvtGen::setup : EvtGen extension not compiled." << std::endl;
	return false;
#endif
}

bool RapidExternalEvtGen::setupGenerator() {
#ifdef RAPID_EVTGEN
	std::cout << "INFO in RapidExternalEvtGen::setupGenerator : Setting up external EvtGen generator." << std::endl;
	EvtRandomEngine* randomEngine = 0;
	EvtAbsRadCorr* radCorrEngine = 0;
	std::list<EvtDecayBase*> extraModels;

	// Define the random number generator
	uint seed = gRandom->GetSeed();
	randomEngine = new EvtMTRandomEngine(seed);

	bool useEvtGenRandom(false);
	EvtExternalGenList genList(true, "", "gamma", useEvtGenRandom);
	radCorrEngine = genList.getPhotosModel();
	extraModels = genList.getListOfModels();

	TString evtPDLPath;
	evtPDLPath += getenv("EVTGEN_ROOT");
	evtPDLPath += "/evt.pdl";

	bool foundDec=false;

	TString decPath;
	decPath += getenv("RAPIDSIM_CONFIG");
	if(decPath!="") {
		decPath += "/config/evtgen/DECAY.DEC";
		if(!gSystem->AccessPathName(decPath)) foundDec=true;
	}

	// We want to initialise EvtGen before we define our DEC file so we can use EvtPDL
	// To do this pass an empty DEC file as the main decay file and pass our file later as a user file
	if(!foundDec) {
		decPath += getenv("RAPIDSIM_ROOT");
		decPath += "/config/evtgen/DECAY.DEC";
	}

	evtGen_ = new EvtGen(decPath.Data(), evtPDLPath.Data(), randomEngine,
			radCorrEngine, &extraModels);

	return true;
#else
	std::cout << "WARNING in RapidExternalEvtGen::setup : EvtGen extension not compiled." << std::endl;
	return false;
#endif
}

void RapidExternalEvtGen::writeDecFile(TString fname, std::vector<RapidParticle*>& parts, bool usePhotos) {
#ifdef RAPID_EVTGEN
	if(!evtGen_) setupGenerator();

	decFileName_ = fname+".DEC";
	std::cout << "INFO in RapidExternalEvtGen::writeDecFile : Writing EvtGen DEC file : " << decFileName_ << std::endl;

	std::ofstream fout;
	fout.open(decFileName_, std::ofstream::out);

	if(usePhotos) {
		fout << "yesPhotos\n" << std::endl;
	} else {
		fout << "noPhotos\n" << std::endl;
	}

	// Loop over all particles and write out Decay rule for each
	for(unsigned int iPart=0; iPart<parts.size(); ++iPart) {
		unsigned int nChildren = parts[iPart]->nDaughters();
		if(nChildren>0) {
			int id = parts[iPart]->id();
			fout << "Decay " << getEvtGenName(id) << "\n1.00\t";
			if ( !(parts[iPart]->evtGenDecayModel()).Contains("TAUOLA") ) {
				for(unsigned int iChild=0; iChild<nChildren; ++iChild) {
					fout << getEvtGenName(parts[iPart]->daughter(iChild)->id()) << "\t";
				}
			}
			fout << parts[iPart]->evtGenDecayModel() << ";" << std::endl;
			fout <<"Enddecay" << std::endl;

			// Workaround to deal with mixing of B0 and Bs
			if(TMath::Abs(id)==531||TMath::Abs(id)==511) fout <<"CDecay " << getEvtGenConjName(id) << std::endl << std::endl;
		}
	}
	fout <<"End\n" << std::endl;
	fout.close();
#else
	std::cout << "WARNING in RapidExternalEvtGen::writeDecFile : EvtGen extension not compiled. Cannot write DEC file "
		  << fname << " for " << parts.size() << "particles with usePhotos=" << usePhotos << "." << std::endl;
#endif
}

TString RapidExternalEvtGen::getEvtGenName(int id) {
#ifdef RAPID_EVTGEN
	EvtId evtId = EvtPDL::evtIdFromLundKC(id);
	TString name = EvtPDL::name(evtId);
	return name;
#else
	std::cout << "WARNING in RapidExternalEvtGen::getEvtGenName : EvtGen extension not compiled. Cannot lookup name for particle ID " << id << "." << std::endl;
	return "";
#endif
}

TString RapidExternalEvtGen::getEvtGenConjName(int id) {
#ifdef RAPID_EVTGEN
	EvtId evtId = EvtPDL::evtIdFromLundKC(id);
	EvtId evtConjId = EvtPDL::chargeConj(evtId);
	TString name = EvtPDL::name(evtConjId);
	return name;
#else
	std::cout << "WARNING in RapidExternalEvtGen::getEvtGenConjName : EvtGen extension not compiled. Cannot lookup conjugate name for particle ID " << id << "." << std::endl;
	return "";
#endif
}
