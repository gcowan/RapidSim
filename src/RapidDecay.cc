#include "RapidDecay.h"

#include "TMath.h"

#include "RooRealVar.h"

#include <fstream>
#include <queue>

#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParticle.h"

RapidDecay::~RapidDecay() {
	if(tree) {
		tree->AutoSave();
		treeFile->Close();
	}
	while(!histos.empty()) {
		delete histos[histos.size()-1];
		histos.pop_back();
	}
	while(!parts.empty()) {
		delete parts[parts.size()-1];
		parts.pop_back();
	}
	std::map<TString, RapidMomentumSmear*>::iterator itr = momSmearCategories.begin();
	while (itr != momSmearCategories.end()) {
		momSmearCategories.erase(itr++);
	}
}

void RapidDecay::loadParentKinematics(TH1F* pt, TH1F* eta) {
	std::cout << "INFO in RapidDecay::loadParentKinematics : setting kinematics of the parent." << std::endl;
	ptHisto=pt;
	etaHisto=eta;
}

bool RapidDecay::loadSmearing(TString category) {
	std::ifstream fin;
	fin.open("../config/smear/"+category, std::ifstream::in);
	if( ! fin.good()) {
		std::cout << "WARNING in RapidDecay::loadSmearing : failed to load smearing category" << category << std::endl;
		fin.close();
		return false;
	}

	TString filename("");
	filename.ReadToken(fin);
	TFile* file = TFile::Open("../rootfiles/smear/"+filename);

	if(!file) {
		std::cout << "WARNING in RapidDecay::loadSmearing : failed to load root file " << filename << std::endl;
		fin.close();
		return false;
	}

	TString type("");
	type.ReadToken(fin);
	if(type=="GAUSS") {
		TString histname("");
		histname.ReadToken(fin);
		TGraphErrors* graph = dynamic_cast<TGraphErrors*>(file->Get(histname));
		if(!graph) {
			std::cout << "WARNING in RapidDecay::loadSmearing : failed to load graph " << histname << std::endl;
			file->Close();
			fin.close();
			return false;
		}

		momSmearCategories[category] = new RapidMomentumSmearGauss(graph);

	} else if(type=="HISTS") {
		double threshold(0.);
		TString histname("");

		std::vector<TH1F*> hists;
		std::vector<double> thresholds;

		while( true ) {
			fin >> threshold;
			histname.ReadToken(fin);
			if(fin.good()) {
				TH1F* hist = dynamic_cast<TH1F*>(file->Get(histname));
				if(hist) {
					hists.push_back(hist);
					thresholds.push_back(threshold);
				} else {
					std::cout << "WARNING in RapidDecay::loadSmearing : failed to load histogram " << histname << std::endl
						  << "                                      threshold will be ignored." << std::endl;

				}
			} else {
				break;
			}
		}

		if(hists.size() == 0) {
			std::cout << "WARNING in RapidDecay::loadSmearing : failed to load any histograms for smearing category " << category << std::endl;
			file->Close();
			fin.close();
			return false;
		}

		momSmearCategories[category] = new RapidMomentumSmearHisto(thresholds, hists);

	} else {
		std::cout << "WARNING in RapidDecay::loadSmearing : unknown smearing type. Category " << category << " not added." << std::endl;
		file->Close();
		fin.close();
		return false;
	}

	fin.close();
	return true;
}

void RapidDecay::setSmearing(unsigned int particle, TString category) {
	if(particle > parts.size()) {
		std::cout << "WARNING in RapidDecay::setSmearing : particle " << particle << " does not exist - smearing functions not set." << std::endl;
		return;
	}
	if(parts[particle]->nDaughters() != 0) {
		std::cout << "WARNING in RapidDecay::setSmearing : particle " << particle << " is composite - smearing functions not set." << std::endl;
		return;
	}

	std::cout << "INFO in RapidDecay::setSmearing : setting smearing functions for particle " << particle << " (category: " << category << ")" << std::endl;

	if(!momSmearCategories.count(category)) {
		if(!loadSmearing(category)) {
			std::cout << "WARNING in RapidDecay::setSmearing : failed to load smearing category " << category << "." << std::endl
				  << "                                     smearing functions not set for particle " << particle << "." << std::endl;
			return;
		}
	}

	parts[particle]->setSmearing(momSmearCategories[category]);
}

void RapidDecay::setAcceptRejectHist(TH1F* hist, ParamType type, std::vector<int> particles) {
	if(accRejHisto) {
		std::cout << "WARNING in RapidDecay::setAcceptRejectHist : accept/reject histogram already set. Original histogram will be used." << std::endl;
		return;
	}

	std::cout << "INFO in RapidDecay::setAcceptRejectHist : setting the required kinematic distribution." << std::endl;
	accRejHisto = hist;
	CustomParameter param;
	param.type = type;
	param.particles = particles;
	accRejParameter = param;

	//correct the histogram to account for the the phasespace distribution
	TH1F* denom = generateAccRejDenominator();
	accRejHisto->Divide(denom);
	delete denom;
}

void RapidDecay::addCustomParameter(TString name, ParamType type, std::vector<int> particles, bool truth, double min, double max) {
	std::cout << "INFO in RapidDecay::addCustomParameter : adding custom parameter " << name << " = ";
	if(truth) std::cout << "TRUE";
	switch(type) {
		case RapidDecay::M:
			std::cout << "M";
			break;
		case RapidDecay::M2:
			std::cout << "M2";
			break;
		case RapidDecay::MT:
			std::cout << "MT";
			break;
		case RapidDecay::E:
			std::cout << "E";
			break;
		case RapidDecay::ET:
			std::cout << "ET";
			break;
		case RapidDecay::P:
			std::cout << "P";
			break;
		case RapidDecay::PX:
			std::cout << "PX";
			break;
		case RapidDecay::PY:
			std::cout << "PY";
			break;
		case RapidDecay::PZ:
			std::cout << "PZ";
			break;
		case RapidDecay::PT:
			std::cout << "PT";
			break;
		case RapidDecay::ETA:
			std::cout << "eta";
			break;
		case RapidDecay::PHI:
			std::cout << "phi";
			break;
		case RapidDecay::RAPIDITY:
			std::cout << "Rapidity";
			break;
		case RapidDecay::GAMMA:
			std::cout << "gamma";
			break;
		case RapidDecay::BETA:
			std::cout << "beta";
			break;
		case RapidDecay::MCORR:
			std::cout << "Mcorr";
			break;
	}
	std::cout << "(";
	for(unsigned int i=0; i<particles.size(); ++i) {
		if(i>0) std::cout << ", ";
		std::cout << parts[particles[i]]->name();
	}
	std::cout << ")" << std::endl;
	CustomParameter param;
	param.name = name;
	param.type = type;
	param.particles = particles;
	param.truth = truth;
	param.minVal = min;
	param.maxVal = max;

	customParams.push_back(param);

	TH1F* hist = new TH1F(name, "", 100, min, max);
	histos.push_back(hist);
}

bool RapidDecay::generate() {
	//keep resonance masses and parent kinematics independent of the accept/reject decision
	//these will only be biased if the function is very inefficient for certain values
	//however, one should not use an a/r function the is highly correlated to these variables
	floatMasses();
	genParent();

	if(accRejHisto) {
		if(!genDecayAccRej()) return false;

	} else {
		if(!genDecay()) return false;
	}

	smearMomenta();
	fillHistos();
	if(tree) fillTree();

	return true;
}

void RapidDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=parts.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		parts[i]->smearMomentum();
	}

}

TLorentzVector RapidDecay::getP(unsigned int i) {
	if(parts.size() > i) return parts[i]->getP();
	else {
		std::cout << "WARNING in RapidDecay::getP : particle: " << i << "does not exist." << std::endl
			  << "                             Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}

TLorentzVector RapidDecay::getPSmeared(unsigned int i) {
	if(parts.size() > i) return parts[i]->getPSmeared();
	else {
		std::cout << "WARNING in RapidDecay::getPSmeared : particle: " << i << "does not exist." << std::endl
			  << "                                    Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}

void RapidDecay::saveHistos(TString fname) {
	std::cout << "INFO in RapidDecay::saveHistos : saving histograms to file: " << fname << std::endl;
	TFile* histFile = new TFile(fname, "RECREATE");

	for(unsigned int i=0; i<histos.size(); ++i) {
		histos[i]->Write();
	}
	histFile->Close();
}

void RapidDecay::saveTree(TString fname) {
	if(!tree) {
		std::cout << "INFO in RapidDecay::saveTree : tree will be saved to file: " << fname << std::endl;
		std::cout << "                            : This will slow down generation." << std::endl;
		treeFile = new TFile(fname, "RECREATE");
		setupTree();
		tree->SetDirectory(treeFile);
	} else {
		std::cout << "WARNING in RapidDecay::saveTree : tree has already been setup and will be saved to the original file." << std::endl;
	}
}

void RapidDecay::loadDecay(TString filename) {
	std::cout << "INFO in RapidDecay::loadDecay : loading decay descriptor from file: " << filename+".decay" << std::endl;
	TString decayStr;
	std::queue<TString> decays;

	std::ifstream fin;
	fin.open(filename+".decay");
	decayStr.ReadLine(fin);
	fin.close();

	decays.push(decayStr);

	while(!decays.empty()) {

		decayStr = decays.front();
		std::cout << "INFO in RapidDecay::loadDecay : Decay descriptor is:" << std::endl
			  << "                                " << decayStr << std::endl;
		decays.pop();

		//first strip out any subdecays and add them to the queue
		while(decayStr.First('{')!=-1) {
			int start = decayStr.Index('{');
			int end = decayStr.Index('}',start);
			if(end < 0) {
				std::cout << "WARNING in RapidDecay::loadDecay : malformed decay descriptor." << std::endl
					  << "                                  Mismatched brackets in:" << decayStr << std::endl;
			}

			//move sub decay into its own string
			TString subDecay = decayStr(start+1,end-start-1);
			subDecay = subDecay.Strip(TString::kBoth);

			//if the sub decay has subdecays of its own then the first '}' won't be the correct one
			while(subDecay.CountChar('{') != subDecay.CountChar('}')) {
				end = decayStr.Index("}",end+1);
				subDecay = decayStr(start+1,end-start-1);
				subDecay = subDecay.Strip(TString::kBoth);
			}
			//add to the queue and label the mother so we know which particle to decay
			decays.push(subDecay);
			std::cout << "INFO in RapidDecay::loadDecay : Found sub-decay:" << std::endl
				  << "                                " << subDecay << std::endl;
			decayStr.Replace(start,end-start+1,"^"+subDecay(0,subDecay.First(" ")));
		}

		//the decaying particle is the first unstable particle with no daughters
		RapidParticle* theMother(0);
		for(unsigned int i=0; i<parts.size(); ++i) {
			RapidParticle* part = parts[i];
			if(!part->stable() && part->nDaughters()==0) {
				theMother = part;
				break;
			}
		}

		//now get info from decay string
		TString token;
		int from(0);

		//first is the parent
		decayStr.Tokenize(token, from, " ");
		//only need to add this if it is the top particle
		if(parts.empty()) {
			theMother = particleData->makeParticle(token, 0);
			parts.push_back(theMother);
		}

		//second should be ->
		decayStr.Tokenize(token, from, " ");
		while (decayStr.Tokenize(token, from, " ")) {
			bool stable = true;
			if(token[0] == '^') {
				token = token.Strip(TString::kBoth,'^');
				stable = false;//flags it to edit later
			}
			RapidParticle* part = particleData->makeParticle(token, theMother);
			part->setStable(stable);
			parts.push_back(part);
			theMother->addDaughter(part);

		}

	}

	loadConfig(filename);

	//automatically add the corrected mass variable if we have any invisible particles
	for(unsigned int i=0; i<parts.size(); ++i) {
		if(parts[i]->invisible()) {
			std::cout << "INFO in RapidDecay::loadDecay : invisible daughter found." << std::endl
				  << "                                Adding corrected mass variable for mother particle." << std::endl;
			std::vector<int> mother;
			mother.push_back(0);
			addCustomParameter("MCorr", RapidDecay::MCORR, mother, false, 0., 7.);
			break;
		}
	}

	setupMasses();
	setupHistos();

	std::cout << "INFO in RapidDecay::loadDecay : Particle summary follows:" << std::endl;
	printf("index\tlabel\t\t   ID\t\tmass (GeV/c^2)\tmother\t\t# daughters\tdaughters\n");
	for(unsigned int i=0; i<parts.size(); ++i) {
		parts[i]->print(i);
	}
}

void RapidDecay::loadConfig(TString filename) {
	std::cout << "INFO in RapidDecay::loadConfig : attempting to load configuration from file: " << filename+".config" << std::endl;

	std::ifstream fin;
	fin.open(filename+".config", std::ifstream::in);
	if( ! fin.good()) {
		std::cout << "INFO in RapidDecay::loadConfig : failed to load configuration. Will write default config file to: " << filename+".config" << std::endl;
		fin.close();
		writeConfig(filename);
		fin.open(filename+".config", std::ifstream::in);
	}

	TString buffer;
	unsigned int currentPart(parts.size());
	while(fin.good()) {
		buffer.ReadLine(fin);
		switch(buffer[0]) {
			case '@': //particle config
				buffer.Remove(TString::kBoth,'@');
				currentPart = buffer.Atoi();
				break;
			case '!': //global config
				break;
			case '#': //comment
				continue;
			default: //continue particle or global config
				int colon = buffer.Index(":");
				TString command = buffer(0,colon);
				command = command.Strip(TString::kBoth);
				TString value = buffer(colon+1, buffer.Length()-colon-1);
				value = value.Strip(TString::kBoth);
				configParticle(currentPart, command, value);
				break;
		}
	}
	std::cout << "INFO in RapidDecay::loadConfig : finished loading configuration." << std::endl;
	fin.close();
}

void RapidDecay::writeConfig(TString filename) {
	std::ofstream fout;
	fout.open(filename+".config", std::ofstream::out);

	for(unsigned int i=0; i<parts.size(); ++i) {
		RapidParticle* part = parts[i];

		fout << "@" << i << "\n";
		fout << "\tname : " << part->name() << "\n";
		if(part->nDaughters()==0) {
			if(TMath::Abs(part->id()) == 11) {
				fout << "\tsmear : electron\n";
			} else {
				fout << "\tsmear : generic\n";
			}
			if(TMath::Abs(part->id()) == 12 ||
					TMath::Abs(part->id()) == 14 ||
					TMath::Abs(part->id()) == 16 ) {
				fout << "\tinvisible : true\n";
			}
		}

	}

	fout.close();
}

void RapidDecay::configParticle(unsigned int part, TString command, TString value) {
	if(part > parts.size()) {
		std::cout << "WARNING in RapidDecay::configParticle : no particle at index " << part << std::endl;
		return;
	}

	if(command=="name") {
		parts[part]->setName(value);
	} else if(command=="smear") {
		setSmearing(part, value);
	} else if(command=="invisible") {
		if(value=="false") {
			parts[part]->setInvisible(false);
		} else if(value=="true") {
			parts[part]->setInvisible(true);
		}
	}
}

void RapidDecay::setupHistos() {
	for(unsigned int i=0; i<parts.size(); ++i) {
		RapidParticle* part = parts[i];

		if(part->nDaughters() < 2) continue;
		// for each subdecay include:
		// parent mass
		// all 2- and 3-body daughter combinations
		// smeared versions of histograms
		TString baseName = part->name();
		TString histName;
		TString axisTitle;
		TH1F* hist(0);

		//for the mother mass
		double mmin = part->minMass();
		double mmax = part->maxMass();
		mmin -= 0.1;
		mmax += 0.1;

		axisTitle = "m(" + baseName + ")";
		histName = baseName+"_M";
		hist = new TH1F(histName, "", 100, mmin, mmax);
		hist->GetXaxis()->SetTitle(axisTitle);
		histos.push_back(hist);

		histName += "_smeared";
		hist = new TH1F(histName, "", 100, mmin, mmax);
		hist->GetXaxis()->SetTitle(axisTitle);
		histos.push_back(hist);

		//2-body IMs
		if(part->nDaughters() > 2) {
			RapidParticle* jDaug = part->daughter(0);

			for( ; jDaug!=0; jDaug=jDaug->next()) {
				for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
					mmin = jDaug->mass() + kDaug->mass();
					mmax = part->mass();
					RapidParticle* daug = part->daughter(0);
					for( ; daug!=0; daug = daug->next()) {
						if(daug==jDaug || daug==kDaug) continue;
						mmax -= daug->mass();
					}
					axisTitle = "m(";
					axisTitle += jDaug->name();
					axisTitle += kDaug->name();
					axisTitle += ")";
					histName = baseName+"_M";
					histName += jDaug->name();
					histName += kDaug->name();
					hist = new TH1F(histName, "", 100, mmin, mmax);
					hist->GetXaxis()->SetTitle(axisTitle);
					histos.push_back(hist);

					histName += "_smeared";
					hist = new TH1F(histName, "", 100, mmin, mmax);
					hist->GetXaxis()->SetTitle(axisTitle);
					histos.push_back(hist);

				}
			}
		}

		//3-body IMs
		if(part->nDaughters() > 3) {
			RapidParticle* jDaug = part->daughter(0);

			for( ; jDaug!=0; jDaug=jDaug->next()) {
				for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
					for(RapidParticle* lDaug=kDaug->next(); lDaug!=0; lDaug=lDaug->next()) {
						mmin = jDaug->mass() + kDaug->mass() + lDaug->mass();
						mmax = part->mass();
						RapidParticle* daug = part->daughter(0);
						for( ; daug!=0; daug = daug->next()) {
							if(daug==jDaug || daug==kDaug || daug==lDaug) continue;
							mmax -= daug->mass();
						}
						axisTitle = "m(";
						axisTitle += jDaug->name();
						axisTitle += kDaug->name();
						axisTitle += lDaug->name();
						axisTitle += ")";
						histName = baseName+"_M";
						histName += jDaug->name();
						histName += kDaug->name();
						histName += lDaug->name();
						hist = new TH1F(histName, "", 100, mmin, mmax);
						hist->GetXaxis()->SetTitle(axisTitle);
						histos.push_back(hist);

						histName += "_smeared";
						hist = new TH1F(histName, "", 100, mmin, mmax);
						hist->GetXaxis()->SetTitle(axisTitle);
						histos.push_back(hist);
					}

				}
			}
		}
	}

}

void RapidDecay::fillHistos() {

	int iHist(0);

	for(unsigned int i=0; i<parts.size(); ++i) {
		RapidParticle* part = parts[i];
		if(part->nDaughters() < 2) continue;
		// for each subdecay include:
		// parent mass
		// all 2- and 3-body daughter combinations
		// smeared versions of histograms

		//for the mother mass
		histos[iHist++]->Fill(getP(i).M());
		histos[iHist++]->Fill(getPSmeared(i).M());

		//2-body IMs
		if(part->nDaughters() > 2) {
			RapidParticle* jDaug = part->daughter(0);

			for( ; jDaug!=0; jDaug=jDaug->next()) {
				for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
					TLorentzVector pSum = jDaug->getP() + kDaug->getP();
					histos[iHist++]->Fill(pSum.M());

					pSum = jDaug->getPSmeared() + kDaug->getPSmeared();
					histos[iHist++]->Fill(pSum.M());
				}
			}
		}

		//3-body IMs
		if(part->nDaughters() > 3) {
			RapidParticle* jDaug = part->daughter(0);

			for( ; jDaug!=0; jDaug=jDaug->next()) {
				for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
					for(RapidParticle* lDaug=kDaug->next(); lDaug!=0; lDaug=lDaug->next()) {
						TLorentzVector pSum = jDaug->getP() + kDaug->getP() + lDaug->getP();
						histos[iHist++]->Fill(pSum.M());

						pSum = jDaug->getPSmeared() + kDaug->getPSmeared() + lDaug->getPSmeared();
						histos[iHist++]->Fill(pSum.M());
					}
				}
			}
		}
	}

	for(unsigned int i=0; i<customParams.size(); ++i) {
		histos[iHist++]->Fill(evalCustomParam(i));
	}
}

void RapidDecay::setupTree() {
	varsPerPart = 19;
	tree = new TTree("DecayTree","DecayTree");
	treeVars = std::vector<double>(parts.size()*varsPerPart + customParams.size(), 0);

	for(unsigned int i=0; i<parts.size(); ++i) {
		TString baseName = parts[i]->name();
		tree->Branch(baseName+"_ID"      ,    &treeVars[varsPerPart*i+0] );
		tree->Branch(baseName+"_M"       ,    &treeVars[varsPerPart*i+1] );
		tree->Branch(baseName+"_E"       ,    &treeVars[varsPerPart*i+2] );
		tree->Branch(baseName+"_P"       ,    &treeVars[varsPerPart*i+3] );
		tree->Branch(baseName+"_PX"      ,    &treeVars[varsPerPart*i+4] );
		tree->Branch(baseName+"_PY"      ,    &treeVars[varsPerPart*i+5] );
		tree->Branch(baseName+"_PZ"      ,    &treeVars[varsPerPart*i+6] );
		tree->Branch(baseName+"_PT"      ,    &treeVars[varsPerPart*i+7] );
		tree->Branch(baseName+"_ETA"     ,    &treeVars[varsPerPart*i+8] );
		tree->Branch(baseName+"_PHI"     ,    &treeVars[varsPerPart*i+9] );
		tree->Branch(baseName+"_TRUEM"   ,    &treeVars[varsPerPart*i+10]);
		tree->Branch(baseName+"_TRUEE"   ,    &treeVars[varsPerPart*i+11]);
		tree->Branch(baseName+"_TRUEP"   ,    &treeVars[varsPerPart*i+12]);
		tree->Branch(baseName+"_TRUEPX"  ,    &treeVars[varsPerPart*i+13]);
		tree->Branch(baseName+"_TRUEPY"  ,    &treeVars[varsPerPart*i+14]);
		tree->Branch(baseName+"_TRUEPZ"  ,    &treeVars[varsPerPart*i+15]);
		tree->Branch(baseName+"_TRUEPT"  ,    &treeVars[varsPerPart*i+16]);
		tree->Branch(baseName+"_TRUEETA" ,    &treeVars[varsPerPart*i+17]);
		tree->Branch(baseName+"_TRUEPHI" ,    &treeVars[varsPerPart*i+18]);
	}
	for(unsigned int i=0; i<customParams.size(); ++i) {
		tree->Branch(customParams[i].name, &treeVars[varsPerPart*parts.size() + i]);
	}
}

void RapidDecay::fillTree() {

	for(unsigned int i=0; i<parts.size(); ++i) {
		TLorentzVector mom = getPSmeared(i);
		treeVars[varsPerPart*i+0] = parts[i]->id();
		treeVars[varsPerPart*i+1] = mom.M();
		treeVars[varsPerPart*i+2] = mom.E();
		treeVars[varsPerPart*i+3] = mom.P();
		treeVars[varsPerPart*i+4] = mom.Px();
		treeVars[varsPerPart*i+5] = mom.Py();
		treeVars[varsPerPart*i+6] = mom.Pz();
		treeVars[varsPerPart*i+7] = mom.Pt();
		treeVars[varsPerPart*i+8] = mom.Eta();
		treeVars[varsPerPart*i+9] = mom.Phi();
		mom = getP(i);
		treeVars[varsPerPart*i+10] = mom.M();
		treeVars[varsPerPart*i+11] = mom.E();
		treeVars[varsPerPart*i+12] = mom.P();
		treeVars[varsPerPart*i+13] = mom.Px();
		treeVars[varsPerPart*i+14] = mom.Py();
		treeVars[varsPerPart*i+15] = mom.Pz();
		treeVars[varsPerPart*i+16] = mom.Pt();
		treeVars[varsPerPart*i+17] = mom.Eta();
		treeVars[varsPerPart*i+18] = mom.Phi();
	}

	for(unsigned int i=0; i<customParams.size(); ++i) {
		treeVars[varsPerPart*parts.size() + i] = evalCustomParam(i);
	}

	tree->Fill();
}

void RapidDecay::floatMasses() {
	for(unsigned int i=0; i<parts.size(); ++i) {
		parts[i]->floatMass();
	}
}

void RapidDecay::setupMasses() {
	for(unsigned int i=0; i<parts.size(); ++i) {
		particleData->setupMass(parts[i]);
	}
}

double RapidDecay::evalCustomParam(int i) {
	CustomParameter param = customParams[i];
	return evalCustomParam(param);
}

double RapidDecay::evalCustomParam(CustomParameter param) {

	if(param.type == RapidDecay::MCORR) return evalCorrectedMass(param);

	TLorentzVector mom;
	if(param.truth) {
		for(unsigned int i=0; i<param.particles.size(); ++i) {
			mom += getP(param.particles[i]);
		}
	} else {
		for(unsigned int i=0; i<param.particles.size(); ++i) {
			mom += getPSmeared(param.particles[i]);
		}
	}
	switch(param.type) {
		case RapidDecay::M:
			return mom.M();
		case RapidDecay::M2:
			return mom.M2();
		case RapidDecay::MT:
			return mom.Mt();
		case RapidDecay::E:
			return mom.E();
		case RapidDecay::ET:
			return mom.Et();
		case RapidDecay::P:
			return mom.P();
		case RapidDecay::PX:
			return mom.Px();
		case RapidDecay::PY:
			return mom.Py();
		case RapidDecay::PZ:
			return mom.Pz();
		case RapidDecay::PT:
			return mom.Pt();
		case RapidDecay::ETA:
			return mom.Eta();
		case RapidDecay::PHI:
			return mom.Phi();
		case RapidDecay::RAPIDITY:
			return mom.Rapidity();
		case RapidDecay::GAMMA:
			return mom.Gamma();
		case RapidDecay::BETA:
			return mom.Beta();
		case RapidDecay::MCORR: //dealt with separately above - included to appease compiler
		default:
			std::cout << "WARNING in RapidDecay::evalCustomParam : unknown parameter type " << param.type << std::endl
				  << "                                         returning 0." << std::endl;

	}

	return 0.;
}

double RapidDecay::evalCorrectedMass(CustomParameter param) {

	TLorentzVector momS, momT;

	//load the true (inc. invisible) and smeared momenta
	for(unsigned int i=0; i<param.particles.size(); ++i) {
		momT += getP(param.particles[i]);
		momS += getPSmeared(param.particles[i]);
	}

	//get the sine and cosine of the angle to the B flight direction with a Gaussian smearing applied
	double cosDir = momT.Vect().Dot(momS.Vect())/(momT.P()*momS.P());
	double dir = TMath::ACos(cosDir) + rand.Gaus(0.,0.01); //TODO smearing parameter should be configurable
	double sinDir = TMath::Sin(dir);
	cosDir = TMath::Cos(dir);

	//get the longitudinal and transverse momentum of the visible daughters wrt the B flight direction
	double pLong  = TMath::Abs(cosDir * momS.P());
	double pTran = TMath::Abs(sinDir * momS.P());

	//invariant masses of the visible daughters and the parent as well as the missing mass
	double mVis2 = momS.M2();
	double mPar2 = momT.M2();
	double mMiss2 = mPar2 - mVis2;

	//the corrected mass
	double mCorr = TMath::Sqrt( mVis2 + pTran*pTran ) + pTran;

	//coefficients of the quadratic equation in pL(invisible)
	double a = 2. * pLong * pLong * mVis2;
	double b = 4*pLong*(2*pTran*pLong - mMiss2);
	double c = 4*pTran*pTran * (pLong*pLong + mPar2) - mMiss2*mMiss2;

	//separate according to whether solutions of pL(invisible) are real or not
	if(b*b - 4.*a*c > 0.) mCorr=-1.*mCorr;

	return mCorr;
}

bool RapidDecay::runAcceptReject() {
	double val = evalCustomParam(accRejParameter);
	double score = accRejHisto->Interpolate(val);
	double max = accRejHisto->GetMaximum();
	if(score > rand.Uniform(max)) return true;
	return false;
}

TH1F* RapidDecay::generateAccRejDenominator() {
	TH1F* denomHisto = dynamic_cast<TH1F*>(accRejHisto->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in RapidDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay(true)) continue;
		denomHisto->Fill(evalCustomParam(accRejParameter));
	}
	return denomHisto;
}

void RapidDecay::genParent() {
	double pt(0), eta(0), phi(rand.Uniform(0,2*TMath::Pi()));
	if(ptHisto)   pt = ptHisto->GetRandom();
	if(etaHisto) eta = etaHisto->GetRandom();
	parts[0]->setPtEtaPhi(pt,eta,phi);
}

bool RapidDecay::genDecay(bool acceptAny) {
	for(unsigned int i=0; i<parts.size(); ++i) {
		RapidParticle* part = parts[i];
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
				while (nGen < maxgen && accept == false){
					accept = event.Generate() > rand.Uniform();
					++nGen;
				} // while

				if(!accept) {
					std::cout << "ERROR in RapidDecay::genDecay : rejected all " << maxgen << " attempts to decay " << part->name() << "." << std::endl;
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

	} while(!passAccRej && ntry<maxgen);

	if(!passAccRej) {
		std::cout << "WARNING in RapidDecay::genDecayAccRej : no events found with required kinematics." << std::endl;
		return false;
	}

	return true;
}
