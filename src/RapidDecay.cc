#include "RapidDecay.h"

#include <fstream>
#include <queue>

#include "TGenPhaseSpace.h"
#include "TMath.h"

#include "RooRealVar.h"

#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParam.h"
#include "RapidParticle.h"

void RapidDecay::loadParentKinematics(TH1F* ptHisto, TH1F* etaHisto) {
	std::cout << "INFO in RapidDecay::loadParentKinematics : setting kinematics of the parent." << std::endl;
	ptHisto_=ptHisto;
	etaHisto_=etaHisto;
}

void RapidDecay::setAcceptRejectHist(TH1F* histo, RapidParam* param) {
	accRejParameter_ = param;
	accRejHisto_     = histo;

	//correct the histogram to account for the the phasespace distribution
	TH1F* denom = generateAccRejDenominator();
	accRejHisto_->Divide(denom);
	delete denom;
}

bool RapidDecay::generate() {
	//keep resonance masses and parent kinematics independent of the accept/reject decision
	//these will only be biased if the function is very inefficient for certain values
	//however, one should not use an a/r function the is highly correlated to these variables
	floatMasses();
	genParent();

	if(accRejHisto_) {
		if(!genDecayAccRej()) return false;

	} else {
		if(!genDecay()) return false;
	}

	smearMomenta();
	writer_.fill();

	return true;
}

void RapidDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=parts_.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		parts_[i]->smearMomentum();
	}

}

void RapidDecay::saveHistos() {
	writer_.save();
}

void RapidDecay::loadDecay(TString filename, bool saveTree) {
	config_.load(filename);

	parts_ = config_.particles();

	TH1F* arHist = config_.acceptRejectHistogram();
	RapidParam* arParam = config_.acceptRejectParameter();

	if(arHist && arParam) {
		setAcceptRejectHist(arHist,arParam);
	}

	setupMasses();
	writer_.setup(parts_, config_.parameters(), filename, saveTree);

	std::cout << "INFO in RapidDecay::loadDecay : Particle summary follows:" << std::endl;
	printf("index\tlabel\t\t   ID\t\tmass (GeV/c^2)\tmother\t\t# daughters\tdaughters\n");
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->print(i);
	}
}

void RapidDecay::floatMasses() {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->floatMass();
	}
}

void RapidDecay::setupMasses() {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		particleData_->setupMass(parts_[i]);
	}
}

bool RapidDecay::runAcceptReject() {
	double val = accRejParameter_->eval();
	int bin = accRejHisto_->FindBin(val);

	double score(0.);
	if(!accRejHisto_->IsBinOverflow(bin) && !accRejHisto_->IsBinUnderflow(bin)) {
		score = accRejHisto_->Interpolate(val);
	}
	double max = accRejHisto_->GetMaximum();
	if(score > gRandom->Uniform(max)) return true;
	return false;
}

TH1F* RapidDecay::generateAccRejDenominator() {
	TH1F* denomHisto = dynamic_cast<TH1F*>(accRejHisto_->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in RapidDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay(true)) continue;
		denomHisto->Fill(accRejParameter_->eval());
	}
	return denomHisto;
}

void RapidDecay::genParent() {
	double pt(0), eta(0), phi(gRandom->Uniform(0,2*TMath::Pi()));
	if(ptHisto_)   pt = ptHisto_->GetRandom();
	if(etaHisto_) eta = etaHisto_->GetRandom();
	parts_[0]->setPtEtaPhi(pt,eta,phi);
}

bool RapidDecay::genDecay(bool acceptAny) {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		if(part->nDaughters()>0) {
			TGenPhaseSpace event;
			// check decay kinematics valid
			if(!event.SetDecay(part->getP(), part->nDaughters(), part->daughterMasses())) {
				std::cout << "ERROR in RapidDecay::genDecay : decay of " << part->name() << " is kinematically forbidden." << std::endl;
				return false;
			}

			// make an event
			if(acceptAny) {
				event.Generate();
			} else {
				int nGen(0);
				bool accept(false);
				while (nGen < maxgen_ && accept == false){
					accept = event.Generate() > gRandom->Uniform();
					++nGen;
				} // while

				if(!accept) {
					std::cout << "ERROR in RapidDecay::genDecay : rejected all " << maxgen_ << " attempts to decay " << part->name() << "." << std::endl;
					return false;
				}
			}

			int j=0;
			for(RapidParticle* jDaug=part->daughter(0); jDaug!=0; jDaug=jDaug->next()) {
				jDaug->setP(*event.GetDecay(j++));
			}
		}
	}

	return true;
}

bool RapidDecay::genDecayAccRej() {
	bool passAccRej(true);
	int ntry(0);

	do {
		if(!genDecay(true)) return false;
		passAccRej = runAcceptReject();
		++ntry;

	} while(!passAccRej && ntry<maxgen_);

	if(!passAccRej) {
		std::cout << "WARNING in RapidDecay::genDecayAccRej : no events found with required kinematics." << std::endl;
		return false;
	}

	return true;
}
