#include "RapidDecay.h"

#include <fstream>
#include <queue>

#include "TGenPhaseSpace.h"
#include "TMath.h"

#include "RooRealVar.h"

#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParticle.h"

RapidDecay::~RapidDecay() {
	if(tree_) {
		tree_->AutoSave();
		treeFile_->Close();
	}
	while(!histos_.empty()) {
		delete histos_[histos_.size()-1];
		histos_.pop_back();
	}
	while(!parts_.empty()) {
		delete parts_[parts_.size()-1];
		parts_.pop_back();
	}
	std::map<TString, RapidMomentumSmear*>::iterator itr = momSmearCategories_.begin();
	while (itr != momSmearCategories_.end()) {
		momSmearCategories_.erase(itr++);
	}
}

void RapidDecay::loadParentKinematics(TH1F* ptHisto, TH1F* etaHisto) {
	std::cout << "INFO in RapidDecay::loadParentKinematics : setting kinematics of the parent." << std::endl;
	ptHisto_=ptHisto;
	etaHisto_=etaHisto;
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

		momSmearCategories_[category] = new RapidMomentumSmearGauss(graph);

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

		momSmearCategories_[category] = new RapidMomentumSmearHisto(thresholds, hists);

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
	if(particle > parts_.size()) {
		std::cout << "WARNING in RapidDecay::setSmearing : particle " << particle << " does not exist - smearing functions not set." << std::endl;
		return;
	}
	if(parts_[particle]->nDaughters() != 0) {
		std::cout << "WARNING in RapidDecay::setSmearing : particle " << particle << " is composite - smearing functions not set." << std::endl;
		return;
	}

	std::cout << "INFO in RapidDecay::setSmearing : setting smearing functions for particle " << particle << " (category: " << category << ")" << std::endl;

	if(!momSmearCategories_.count(category)) {
		if(!loadSmearing(category)) {
			std::cout << "WARNING in RapidDecay::setSmearing : failed to load smearing category " << category << "." << std::endl
				  << "                                     smearing functions not set for particle " << particle << "." << std::endl;
			return;
		}
	}

	parts_[particle]->setSmearing(momSmearCategories_[category]);
}

void RapidDecay::setAcceptRejectHist(TString histFile, TString histName, RapidParam* param) {
	if(accRejHisto_) {
		std::cout << "WARNING in RapidDecay::setAcceptRejectHist : accept/reject histogram already set." << std::endl
			  << "                                             original histogram will be used." << std::endl;
		return;
	}

	TFile* file = TFile::Open(histFile);
	if(!file) {
		std::cout << "WARNING in RapidDecay::setAcceptRejectHist : could not open file " << histFile << "." << std::endl
			  << "                                             accept/reject histogram not set." << std::endl;
		return;
	}
	accRejHisto_ = dynamic_cast<TH1F*>(file->Get(histName));
	if(!accRejHisto_) {
		std::cout << histName << std::endl;
		std::cout << "WARNING in RapidDecay::setAcceptRejectHist : could not load histogram " << histName << "." << std::endl
			  << "                                             accept/reject histogram not set." << std::endl;
		return;
	}

	std::cout << "INFO in RapidDecay::setAcceptRejectHist : setting the required kinematic distribution." << std::endl;
	accRejParameter_ = param;

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
	fillHistos();
	if(tree_) fillTree();

	return true;
}

void RapidDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=parts_.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		parts_[i]->smearMomentum();
	}

}

TLorentzVector RapidDecay::getP(unsigned int i) {
	if(parts_.size() > i) return parts_[i]->getP();
	else {
		std::cout << "WARNING in RapidDecay::getP : particle: " << i << "does not exist." << std::endl
			  << "                             Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}

TLorentzVector RapidDecay::getPSmeared(unsigned int i) {
	if(parts_.size() > i) return parts_[i]->getPSmeared();
	else {
		std::cout << "WARNING in RapidDecay::getPSmeared : particle: " << i << "does not exist." << std::endl
			  << "                                    Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}

void RapidDecay::saveHistos(TString fname) {
	std::cout << "INFO in RapidDecay::saveHistos : saving histograms to file: " << fname << std::endl;
	TFile* histFile = new TFile(fname, "RECREATE");

	for(unsigned int i=0; i<histos_.size(); ++i) {
		histos_[i]->Write();
	}
	histFile->Close();
}

void RapidDecay::saveTree(TString fname) {
	if(!tree_) {
		std::cout << "INFO in RapidDecay::saveTree : tree will be saved to file: " << fname << std::endl;
		std::cout << "                            : This will slow down generation." << std::endl;
		treeFile_ = new TFile(fname, "RECREATE");
		setupTree();
		tree_->SetDirectory(treeFile_);
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
		for(unsigned int i=0; i<parts_.size(); ++i) {
			RapidParticle* part = parts_[i];
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
		if(parts_.empty()) {
			theMother = particleData_->makeParticle(token, 0);
			parts_.push_back(theMother);
		}

		//second should be ->
		decayStr.Tokenize(token, from, " ");
		while (decayStr.Tokenize(token, from, " ")) {
			bool stable = true;
			if(token[0] == '^') {
				token = token.Strip(TString::kBoth,'^');
				stable = false;//flags it to edit later
			}
			RapidParticle* part = particleData_->makeParticle(token, theMother);
			part->setStable(stable);
			parts_.push_back(part);
			theMother->addDaughter(part);

		}

	}

	loadConfig(filename);

	//automatically add the corrected mass variable if we have any invisible particles
	for(unsigned int i=0; i<parts_.size(); ++i) {
		if(parts_[i]->invisible()) {
			std::cout << "INFO in RapidDecay::loadDecay : invisible daughter found." << std::endl
				  << "                                Adding corrected mass variable for mother particle." << std::endl;
			std::vector<int> mother;
			mother.push_back(0);
			RapidParam* param = new RapidParam("MCorr", RapidParam::MCORR, parts_[0], false);
			customParams_.push_back(param);
			break;
		}
	}

	setupMasses();
	setupHistos();

	std::cout << "INFO in RapidDecay::loadDecay : Particle summary follows:" << std::endl;
	printf("index\tlabel\t\t   ID\t\tmass (GeV/c^2)\tmother\t\t# daughters\tdaughters\n");
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->print(i);
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
	unsigned int currentPart(parts_.size());
	while(fin.good()) {
		buffer.ReadLine(fin);
		switch(buffer[0]) {
			case '@': //particle config
				buffer.Remove(TString::kBoth,'@');
				currentPart = buffer.Atoi();
				break;
			case '!': //global config
				currentPart = parts_.size();
				break;
			case '#': //comment
				continue;
			default: //continue particle or global config
				int colon = buffer.Index(":");
				TString command = buffer(0,colon);
				command = command.Strip(TString::kBoth);
				TString value = buffer(colon+1, buffer.Length()-colon-1);
				value = value.Strip(TString::kBoth);

				if(currentPart==parts_.size()) configGlobal(command, value);
				else configParticle(currentPart, command, value);
				break;
		}
	}
	std::cout << "INFO in RapidDecay::loadConfig : finished loading configuration." << std::endl;
	fin.close();
}

void RapidDecay::writeConfig(TString filename) {
	std::ofstream fout;
	fout.open(filename+".config", std::ofstream::out);

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];

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
	if(part > parts_.size()) {
		std::cout << "WARNING in RapidDecay::configParticle : no particle at index " << part << std::endl;
		return;
	}

	if(command=="name") {
		parts_[part]->setName(value);
	} else if(command=="smear") {
		setSmearing(part, value);
	} else if(command=="invisible") {
		if(value=="false") {
			parts_[part]->setInvisible(false);
		} else if(value=="true") {
			parts_[part]->setInvisible(true);
		}
	}
}

void RapidDecay::configGlobal(TString command, TString value) {
	if(command=="seed") {
		int seed = value.Atoi();
		gRandom->SetSeed(seed);
		std::cout << "INFO in RapidDecay::configGlobal : setting seed for random number generation to " << seed << "." << std::endl
		          << "                                   seed is now " << gRandom->GetSeed() << "." << std::endl;
	} else if(command=="param") {
		RapidParam* param = loadParam(value);
		if(param) customParams_.push_back(param);
	} else if(command=="shape") {
		int from(0);
		TString histName, histFile;
		value.Tokenize(histFile,from," ");
		value.Tokenize(histName,from," ");
		TString paramStr = value(from, value.Length()-from);
		histFile = histFile.Strip(TString::kBoth);
		histName = histName.Strip(TString::kBoth);
		paramStr = paramStr.Strip(TString::kBoth);

		RapidParam* param = loadParam(paramStr);
		if(param) setAcceptRejectHist(histFile, histName, param);
	}
}

RapidParam* RapidDecay::loadParam(TString paramStr) {
	int from(0);
	TString name, buffer;
	RapidParam::ParamType type;
	bool truth = false;
	std::vector<RapidParticle*> partlist;
	paramStr.Tokenize(name,from," ");
	paramStr.Tokenize(buffer,from," ");
	type = RapidParam::typeFromString(buffer);

	while(paramStr.Tokenize(buffer,from," ")) {
		if(buffer.IsDec()) {
			unsigned int index = buffer.Atoi();
			if(index>=parts_.size()) {
				std::cout << "WARNING in RapidDecay::loadParam : particle " << index << "is out of range." << std::endl
					  << "                                   parameter " << name << " has not been added." << std::endl;
				return 0;
			} else {
				partlist.push_back(parts_[index]);
			}
		} else if(buffer=="TRUE") {
			truth = true;
		} else {
			std::cout << "WARNING in RapidDecay::loadParam : unknown argument " << buffer << "in parameter " << name << "." << std::endl
				  << "                                   argument will be ignored." << std::endl;
		}
	}

	if(type == RapidParam::THETA && partlist.size()!=2) {
		std::cout << "WARNING in RapidDecay::loadParam : theta parameter expects two particles, " << partlist.size() << "were given." << std::endl
			  << "                                   parameter " << name << " has not been added." << std::endl;
		return 0;
	}

	RapidParam* param = new RapidParam(name,type,partlist,truth);
	return param;
}

void RapidDecay::setupHistos() {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];

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
		histos_.push_back(hist);

		histName += "_smeared";
		hist = new TH1F(histName, "", 100, mmin, mmax);
		hist->GetXaxis()->SetTitle(axisTitle);
		histos_.push_back(hist);

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
					histos_.push_back(hist);

					histName += "_smeared";
					hist = new TH1F(histName, "", 100, mmin, mmax);
					hist->GetXaxis()->SetTitle(axisTitle);
					histos_.push_back(hist);

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
						histos_.push_back(hist);

						histName += "_smeared";
						hist = new TH1F(histName, "", 100, mmin, mmax);
						hist->GetXaxis()->SetTitle(axisTitle);
						histos_.push_back(hist);
					}

				}
			}
		}
	}

	std::vector<RapidParam*>::iterator it = customParams_.begin();
	for( ; it!= customParams_.end(); ++it) {
		RapidParam* param = (*it);
		TH1F* hist = new TH1F(param->name(), "", 100, param->min(), param->max());
		histos_.push_back(hist);
	}

}

void RapidDecay::fillHistos() {

	int iHist(0);

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		if(part->nDaughters() < 2) continue;
		// for each subdecay include:
		// parent mass
		// all 2- and 3-body daughter combinations
		// smeared versions of histograms

		//for the mother mass
		histos_[iHist++]->Fill(getP(i).M());
		histos_[iHist++]->Fill(getPSmeared(i).M());

		//2-body IMs
		if(part->nDaughters() > 2) {
			RapidParticle* jDaug = part->daughter(0);

			for( ; jDaug!=0; jDaug=jDaug->next()) {
				for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
					TLorentzVector pSum = jDaug->getP() + kDaug->getP();
					histos_[iHist++]->Fill(pSum.M());

					pSum = jDaug->getPSmeared() + kDaug->getPSmeared();
					histos_[iHist++]->Fill(pSum.M());
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
						histos_[iHist++]->Fill(pSum.M());

						pSum = jDaug->getPSmeared() + kDaug->getPSmeared() + lDaug->getPSmeared();
						histos_[iHist++]->Fill(pSum.M());
					}
				}
			}
		}
	}

	std::vector<RapidParam*>::iterator it = customParams_.begin();
	for( ; it!= customParams_.end(); ++it) {
		histos_[iHist++]->Fill((*it)->eval());
	}
}

void RapidDecay::setupTree() {
	varsPerPart_ = 19;
	tree_ = new TTree("DecayTree","DecayTree");
	treeVars_ = std::vector<double>(parts_.size()*varsPerPart_ + customParams_.size(), 0);

	for(unsigned int i=0; i<parts_.size(); ++i) {
		TString baseName = parts_[i]->name();
		tree_->Branch(baseName+"_ID"      ,    &treeVars_[varsPerPart_*i+0] );
		tree_->Branch(baseName+"_M"       ,    &treeVars_[varsPerPart_*i+1] );
		tree_->Branch(baseName+"_E"       ,    &treeVars_[varsPerPart_*i+2] );
		tree_->Branch(baseName+"_P"       ,    &treeVars_[varsPerPart_*i+3] );
		tree_->Branch(baseName+"_PX"      ,    &treeVars_[varsPerPart_*i+4] );
		tree_->Branch(baseName+"_PY"      ,    &treeVars_[varsPerPart_*i+5] );
		tree_->Branch(baseName+"_PZ"      ,    &treeVars_[varsPerPart_*i+6] );
		tree_->Branch(baseName+"_PT"      ,    &treeVars_[varsPerPart_*i+7] );
		tree_->Branch(baseName+"_ETA"     ,    &treeVars_[varsPerPart_*i+8] );
		tree_->Branch(baseName+"_PHI"     ,    &treeVars_[varsPerPart_*i+9] );
		tree_->Branch(baseName+"_TRUEM"   ,    &treeVars_[varsPerPart_*i+10]);
		tree_->Branch(baseName+"_TRUEE"   ,    &treeVars_[varsPerPart_*i+11]);
		tree_->Branch(baseName+"_TRUEP"   ,    &treeVars_[varsPerPart_*i+12]);
		tree_->Branch(baseName+"_TRUEPX"  ,    &treeVars_[varsPerPart_*i+13]);
		tree_->Branch(baseName+"_TRUEPY"  ,    &treeVars_[varsPerPart_*i+14]);
		tree_->Branch(baseName+"_TRUEPZ"  ,    &treeVars_[varsPerPart_*i+15]);
		tree_->Branch(baseName+"_TRUEPT"  ,    &treeVars_[varsPerPart_*i+16]);
		tree_->Branch(baseName+"_TRUEETA" ,    &treeVars_[varsPerPart_*i+17]);
		tree_->Branch(baseName+"_TRUEPHI" ,    &treeVars_[varsPerPart_*i+18]);
	}
	for(unsigned int i=0; i<customParams_.size(); ++i) {
		tree_->Branch(customParams_[i]->name(), &treeVars_[varsPerPart_*parts_.size() + i]);
	}
}

void RapidDecay::fillTree() {

	for(unsigned int i=0; i<parts_.size(); ++i) {
		TLorentzVector mom = getPSmeared(i);
		treeVars_[varsPerPart_*i+0] = parts_[i]->id();
		treeVars_[varsPerPart_*i+1] = mom.M();
		treeVars_[varsPerPart_*i+2] = mom.E();
		treeVars_[varsPerPart_*i+3] = mom.P();
		treeVars_[varsPerPart_*i+4] = mom.Px();
		treeVars_[varsPerPart_*i+5] = mom.Py();
		treeVars_[varsPerPart_*i+6] = mom.Pz();
		treeVars_[varsPerPart_*i+7] = mom.Pt();
		treeVars_[varsPerPart_*i+8] = mom.Eta();
		treeVars_[varsPerPart_*i+9] = mom.Phi();
		mom = getP(i);
		treeVars_[varsPerPart_*i+10] = mom.M();
		treeVars_[varsPerPart_*i+11] = mom.E();
		treeVars_[varsPerPart_*i+12] = mom.P();
		treeVars_[varsPerPart_*i+13] = mom.Px();
		treeVars_[varsPerPart_*i+14] = mom.Py();
		treeVars_[varsPerPart_*i+15] = mom.Pz();
		treeVars_[varsPerPart_*i+16] = mom.Pt();
		treeVars_[varsPerPart_*i+17] = mom.Eta();
		treeVars_[varsPerPart_*i+18] = mom.Phi();
	}

	for(unsigned int i=0; i<customParams_.size(); ++i) {
		treeVars_[varsPerPart_*parts_.size() + i] = customParams_[i]->eval();
	}

	tree_->Fill();
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
