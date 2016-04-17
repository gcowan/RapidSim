#include "RapidConfig.h"

#include <fstream>
#include <iostream>
#include <queue>

#include "TFile.h"
#include "TRandom.h"

#include "RapidAcceptance.h"
#include "RapidDecay.h"
#include "RapidHistWriter.h"
#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParam.h"
#include "RapidParticle.h"
#include "RapidParticleData.h"

RapidConfig::~RapidConfig() {
	std::map<TString, RapidMomentumSmear*>::iterator itr = momSmearCategories_.begin();
	while (itr != momSmearCategories_.end()) {
		momSmearCategories_.erase(itr++);
	}

	while(!parts_.empty()) {
		delete parts_[parts_.size()-1];
		parts_.pop_back();
	}

	while(!params_.empty()) {
		delete params_[params_.size()-1];
		params_.pop_back();
	}

	if(accRejHisto_) delete accRejHisto_;
	if(accRejParameter_) delete accRejParameter_;
}

bool RapidConfig::load(TString fileName) {
	fileName_ = fileName;

	if(!loadDecay()) return false;
	if(!loadConfig()) return false;

	//automatically add the corrected mass variable if we have any invisible particles
	for(unsigned int i=0; i<parts_.size(); ++i) {
		if(parts_[i]->invisible()) {
			std::cout << "INFO in RapidConfig::load : invisible daughter found." << std::endl
				  << "                            Adding corrected mass variable for mother particle." << std::endl;
			std::vector<int> mother;
			mother.push_back(0);
			RapidParam* param = new RapidParam("MCorr", RapidParam::MCORR, parts_[0], false);
			params_.push_back(param);
			break;
		}
	}

	return true;
}

RapidDecay* RapidConfig::getDecay() {
	if(!decay_) {
		//check we have a particle to decay
		if(parts_.empty()) return 0;

		decay_ = new RapidDecay(parts_);

		//check that we have an FONLL model for the parent kinematics
		if(!loadParentKinematics()) {
			return 0;
		}
		decay_->setParentKinematics(ptHisto_,etaHisto_);

		//load any PDF to generate
		if(accRejHisto_ && accRejParameter_) {
			decay_->setAcceptRejectHist(accRejHisto_,accRejParameter_);
		}
	}

	return decay_;
}

RapidAcceptance* RapidConfig::getAcceptance() {
	if(!acceptance_) {
		acceptance_ = new RapidAcceptance(acceptanceType_, parts_);
	}
	return acceptance_;
}

RapidHistWriter* RapidConfig::getWriter(bool saveTree) {
	if(!writer_) {
		writer_ = new RapidHistWriter(parts_, params_, fileName_, saveTree);
	}

	return writer_;
}

bool RapidConfig::loadDecay() {
	std::cout << "INFO in RapidConfig::loadDecay : loading decay descriptor from file: " << fileName_ << ".decay" << std::endl;
	TString decayStr;
	std::queue<TString> decays;
	std::queue<RapidParticle*> mothers;

	std::ifstream fin;
	fin.open(fileName_+".decay");
	if(!fin.good()) {
		std::cout << "ERROR in RapidConfig::loadDecay : file " << fileName_ << ".decay not found." << std::endl;
		return false;
	}
	decayStr.ReadLine(fin);
	fin.close();

	std::cout << "INFO in RapidConfig::loadDecay : Decay descriptor is:" << std::endl
		  << "                                 " << decayStr << std::endl;

	decays.push(decayStr);

	RapidParticleData* particleData = RapidParticleData::getInstance();

	while(!decays.empty()) {

		decayStr = decays.front();
		decays.pop();

		//first strip out any subdecays and add them to the queue
		while(decayStr.First('{')!=-1) {
			int start = decayStr.Index('{');
			int end = decayStr.Index('}',start);
			if(end < 0) {
				std::cout << "ERROR in RapidConfig::loadDecay : malformed decay descriptor." << std::endl
					  << "                                  Mismatched brackets in:" << decayStr << std::endl;
				return false;
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
			std::cout << "INFO in RapidConfig::loadDecay : Found sub-decay:" << std::endl
				  << "                                 " << subDecay << std::endl;
			decayStr.Replace(start,end-start+1,"^"+subDecay(0,subDecay.First(" ")));
		}

		//the decaying particle is the first particle remaining in the mothers list
		RapidParticle* theMother(0);
		if(!mothers.empty()) {
			theMother = mothers.front();
			mothers.pop();
		}

		//now get info from decay string
		TString token;
		int from(0);

		//first is the parent
		decayStr.Tokenize(token, from, " ");
		//only need to add this if it is the top particle
		if(parts_.empty()) {
			theMother = particleData->makeParticle(token, 0);
			parts_.push_back(theMother);

			//get flavour of mother for FONLL
			if(theMother->hasBeauty()) {
				motherFlavour_ = "b";
			} else if(theMother->hasCharm()) {
				motherFlavour_ = "c";
			} else {
				std::cout << "WARNING in RapidConfig::loadDecay : Mother has neither beauty nor charm." << std::endl;
				std::cout << "                                    defaulting to b-quark kinematics." << std::endl;
				motherFlavour_ = "b";
			}
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
			if(!stable) mothers.push(part);
			parts_.push_back(part);
			theMother->addDaughter(part);

		}

	}

	return true;
}

bool RapidConfig::loadConfig() {
	std::cout << "INFO in RapidConfig::loadConfig : attempting to load configuration from file: " << fileName_+".config" << std::endl;

	std::ifstream fin;
	fin.open(fileName_+".config", std::ifstream::in);
	if( ! fin.good()) {
		std::cout << "INFO in RapidConfig::loadConfig : failed to load configuration. Will write default config file to: " << fileName_+".config" << std::endl;
		fin.close();
		writeConfig();
		fin.open(fileName_+".config", std::ifstream::in);
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
	std::cout << "INFO in RapidConfig::loadConfig : finished loading configuration." << std::endl;
	fin.close();

	return true;
}

void RapidConfig::writeConfig() {
	std::ofstream fout;
	fout.open(fileName_+".config", std::ofstream::out);

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

void RapidConfig::configParticle(unsigned int part, TString command, TString value) {
	if(part > parts_.size()) {
		std::cout << "WARNING in RapidConfig::configParticle : no particle at index " << part << std::endl;
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

void RapidConfig::configGlobal(TString command, TString value) {
	if(command=="seed") {
		int seed = value.Atoi();
		gRandom->SetSeed(seed);
		std::cout << "INFO in RapidConfig::configGlobal : setting seed for random number generation to " << seed << "." << std::endl
		          << "                                    seed is now " << gRandom->GetSeed() << "." << std::endl;
	} else if(command=="acceptance") {
		acceptanceType_ = RapidAcceptance::typeFromString(value);
	} else if(command=="energy") {
		ppEnergy_ = value.Atof();
	} else if(command=="mother") {
		motherFlavour_ = value;
	} else if(command=="param") {
		RapidParam* param = loadParam(value);
		if(param) params_.push_back(param);
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
		if(param) loadAcceptRejectHist(histFile, histName, param);
	}
}

RapidParam* RapidConfig::loadParam(TString paramStr) {
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
				std::cout << "WARNING in RapidConfig::loadParam : particle " << index << "is out of range." << std::endl
					  << "                                    parameter " << name << " has not been added." << std::endl;
				return 0;
			} else {
				partlist.push_back(parts_[index]);
			}
		} else if(buffer=="TRUE") {
			truth = true;
		} else {
			std::cout << "WARNING in RapidConfig::loadParam : unknown argument " << buffer << "in parameter " << name << "." << std::endl
				  << "                                    argument will be ignored." << std::endl;
		}
	}

	if(type == RapidParam::THETA && partlist.size()!=2) {
		std::cout << "WARNING in RapidConfig::loadParam : theta parameter expects two particles, " << partlist.size() << "were given." << std::endl
			  << "                                    parameter " << name << " has not been added." << std::endl;
		return 0;
	}

	RapidParam* param = new RapidParam(name,type,partlist,truth);
	return param;
}

bool RapidConfig::loadSmearing(TString category) {
	std::ifstream fin;
	fin.open("../config/smear/"+category, std::ifstream::in);
	if( ! fin.good()) {
		std::cout << "WARNING in RapidConfig::loadSmearing : failed to load smearing category" << category << std::endl;
		fin.close();
		return false;
	}

	TString filename("");
	filename.ReadToken(fin);
	TFile* file = TFile::Open("../rootfiles/smear/"+filename);

	if(!file) {
		std::cout << "WARNING in RapidConfig::loadSmearing : failed to load root file " << filename << std::endl;
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
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load graph " << histname << std::endl;
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
					std::cout << "WARNING in RapidConfig::loadSmearing : failed to load histogram " << histname << std::endl
						  << "                                       threshold will be ignored." << std::endl;

				}
			} else {
				break;
			}
		}

		if(hists.size() == 0) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load any histograms for smearing category " << category << std::endl;
			file->Close();
			fin.close();
			return false;
		}

		momSmearCategories_[category] = new RapidMomentumSmearHisto(thresholds, hists);

	} else {
		std::cout << "WARNING in RapidConfig::loadSmearing : unknown smearing type. Category " << category << " not added." << std::endl;
		file->Close();
		fin.close();
		return false;
	}

	fin.close();
	return true;
}

void RapidConfig::setSmearing(unsigned int particle, TString category) {
	if(particle > parts_.size()) {
		std::cout << "WARNING in RapidConfig::setSmearing : particle " << particle << " does not exist - smearing functions not set." << std::endl;
		return;
	}
	if(parts_[particle]->nDaughters() != 0) {
		std::cout << "WARNING in RapidConfig::setSmearing : particle " << particle << " is composite - smearing functions not set." << std::endl;
		return;
	}

	std::cout << "INFO in RapidConfig::setSmearing : setting smearing functions for particle " << particle << " (category: " << category << ")" << std::endl;

	if(!momSmearCategories_.count(category)) {
		if(!loadSmearing(category)) {
			std::cout << "WARNING in RapidConfig::setSmearing : failed to load smearing category " << category << "." << std::endl
				  << "                                      smearing functions not set for particle " << particle << "." << std::endl;
			return;
		}
	}

	parts_[particle]->setSmearing(momSmearCategories_[category]);
}

void RapidConfig::loadAcceptRejectHist(TString histFile, TString histName, RapidParam* param) {
	if(accRejHisto_) {
		std::cout << "WARNING in RapidConfig::loadAcceptRejectHist : accept/reject histogram already set." << std::endl
			  << "                                               original histogram will be used." << std::endl;
		return;
	}

	TFile* file = TFile::Open(histFile);
	if(!file) {
		std::cout << "WARNING in RapidConfig::loadAcceptRejectHist : could not open file " << histFile << "." << std::endl
			  << "                                               accept/reject histogram not set." << std::endl;
		return;
	}
	accRejHisto_ = dynamic_cast<TH1F*>(file->Get(histName));
	if(!accRejHisto_) {
		std::cout << "WARNING in RapidConfig::loadAcceptRejectHist : could not load histogram " << histName << "." << std::endl
			  << "                                               accept/reject histogram not set." << std::endl;
		return;
	}

	std::cout << "INFO in RapidConfig::loadAcceptRejectHist : setting the required kinematic distribution." << std::endl;
	accRejParameter_ = param;
}

bool RapidConfig::loadParentKinematics() {
	TString fileName("../rootfiles/fonll/lhc_");
	fileName += motherFlavour_; fileName += "_";
	fileName += ppEnergy_; fileName += "tev";
	fileName += ".root";
	TFile* file = TFile::Open(fileName);

	if(!file) {
		std::cout << "ERROR in RapidConfig::loadParentKinematics : unknown kinematics " << motherFlavour_ << "-quark from " << ppEnergy_ << " TeV pp collision." << std::endl;
		std::cout << "                                             file " << fileName << " not found." << std::endl;
		return false;
	}

	ptHisto_ = dynamic_cast<TH1F*>(file->Get("pthisto"));
	etaHisto_ = dynamic_cast<TH1F*>(file->Get("etahisto"));

	return true;
}
