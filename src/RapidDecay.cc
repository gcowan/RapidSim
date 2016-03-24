#include "RapidDecay.h"

#include "TMath.h"

#include "RooRealVar.h"

#include <fstream>
#include <queue>

#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"

RapidDecay::~RapidDecay() {
	if(tree) {
		tree->AutoSave();
		treeFile->Close();
	}
	while(!histos.empty()) {
		delete histos[histos.size()-1];
		histos.pop_back();
	}
	std::map<int, RapidMomentumSmear*>::iterator itr = momSmear.begin();
	while (itr != momSmear.end()) {
		momSmear.erase(itr++);
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
	if(nDaug[particle] != 0) {
		std::cout << "WARNING in RapidDecay::setSmearing : particle " << particle << " is composite - smearing functions not set." << std::endl;
		return;
	}
	if(momSmear.count(particle)) {
		std::cout << "WARNING in RapidDecay::setSmearing : smearing function for particle " << particle << " has already been set." << std::endl;
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

	momSmear[particle] = momSmearCategories[category];
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
		std::cout << names[particles[i]];
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

	bool passAccRej(true);
	int ntry(0);

	do {
		if(!genDecay()) return false;
		if(accRejHisto) passAccRej = runAcceptReject();
		++ntry;

	} while(!passAccRej && ntry<maxgen);

	if(!passAccRej) {
		std::cout << "WARNING in RapidDecay::generate : no events found with required kinematics." << std::endl;
		return false;
	}

	smearMomenta();
	fillHistos();
	if(tree) fillTree();

	return true;
}

void RapidDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=p.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		if(nDaug[i] == 0) {
			//smear momentum of each detected particle
			if(invisible[i]) {
				pSmeared[i].SetXYZM(0.,0.,0.,0.);
			} else if(momSmear.count(i)) {
				pSmeared[i] = momSmear[i]->smearMomentum(p[i]);
			} else {
				pSmeared[i] = p[i];
			}
			//pSmeared[i] = smearedVec(p[i],smearGraph,rand);
		} else {
			//reconstruct mothers from their daughters
			pSmeared[i] = TLorentzVector();
			int first = firstDaug[i];
			for(int j=0; j<nDaug[i]; ++j) {
				pSmeared[i] += pSmeared[first+j];
			}
		}
	}

}

TLorentzVector RapidDecay::getP(unsigned int i) {
	if(p.size() > i) return p[i];
	else {
		std::cout << "WARNING in RapidDecay::getP : particle: " << i << "does not exist." << std::endl
			<< "                             Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}

TLorentzVector RapidDecay::getPSmeared(unsigned int i) {
	if(pSmeared.size() > i) return pSmeared[i];
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

		//the decaying particle is the first one in parts with nDaug set to -1 
		int theMother(0);
		for(unsigned int i=0; i<nDaug.size(); ++i) {
			if(nDaug[i]==-1) {
				theMother = i;
				break;
			}
		}

		//now get info from decay string
		TString token;
		int from(0);

		int count(0);

		//first is the parent
		decayStr.Tokenize(token, from, " ");
		if(parts.empty()) {
			int id = pdgCode(token);
			names.push_back(getUniqName(token));
			parts.push_back(id);
			mother.push_back(-1);
			nDaug.push_back(-1);//-1 flags it to edit later
			firstDaug.push_back(1);
			m.push_back(getMass(id));
			p.push_back(TLorentzVector());
			pSmeared.push_back(TLorentzVector());
			invisible.push_back(false);
		}

		//second should be ->
		decayStr.Tokenize(token, from, " ");
		while (decayStr.Tokenize(token, from, " ")) {
			if(token[0] == '^') {
				token = token.Strip(TString::kBoth,'^');
				nDaug.push_back(-1);//-1 flags it to edit later
			} else {
				nDaug.push_back(0);
			}
			int id = pdgCode(token);
			names.push_back(getUniqName(token));
			parts.push_back(id);
			mother.push_back(theMother);
			firstDaug.push_back(-1);
			m.push_back(getMass(id));
			p.push_back(TLorentzVector());
			pSmeared.push_back(TLorentzVector());
			invisible.push_back(false);
			++count;

		}

		nDaug[theMother] = count;
		firstDaug[theMother] = parts.size() - count;

	}

	loadConfig(filename);

	//automatically add the corrected mass variable if we have any invisible particles
	for(unsigned int i=0; i<parts.size(); ++i) {
		if(invisible[i]) {
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
	printf("index\tlabel\t\t   ID\t\tmass (GeV/c^2)\tmother\t# daughters\tfirst daughter\n");
	for(unsigned int i=0; i<parts.size(); ++i) {
		printf("%3d\t%-15s\t%6d\t\t%.6f\t%3d\t%2d\t\t%3d\n", i, names[i].Data(), parts[i], m[i], mother[i], nDaug[i], firstDaug[i]);
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
		fout << "@" << i << "\n";
		fout << "\tname : " << names[i] << "\n";
		if(nDaug[i]==0) {
			if(TMath::Abs(parts[i]) == 11) {
				fout << "\tsmear : electron\n";
			} else {
				fout << "\tsmear : generic\n";
			}
			if(TMath::Abs(parts[i]) == 12 ||
			   TMath::Abs(parts[i]) == 14 ||
			   TMath::Abs(parts[i]) == 16 ) {
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
		names[part] = value;
	} else if(command=="smear") {
		setSmearing(part, value);
	} else if(command=="invisible") {
		if(value=="false") {
			invisible[part] = false;
		} else if(value=="true") {
			invisible[part] = true;
		}
	}
}

TString RapidDecay::getUniqName(TString base) {
	//first sanitise
	base = base.ReplaceAll("+","p");
	base = base.ReplaceAll("-","m");
	base = base.ReplaceAll("*","st");
	base = base.ReplaceAll("(","_");
	base = base.ReplaceAll(")","_");
	base = base.ReplaceAll("[","_");
	base = base.ReplaceAll("]","_");
	base = base.ReplaceAll("<","_");
	base = base.ReplaceAll(">","_");
	base = base.ReplaceAll("{","_");
	base = base.ReplaceAll("}","_");
	base = base.ReplaceAll(" ","_");
	base = base.ReplaceAll("$","");
	base = base.ReplaceAll("%","");
	base = base.ReplaceAll("&","");
	base = base.ReplaceAll("/","");
	base = base.ReplaceAll(":","");
	base = base.ReplaceAll(";","");
	base = base.ReplaceAll("=","");
	base = base.ReplaceAll("\\","");
	base = base.ReplaceAll("^","");
	base = base.ReplaceAll("|","");
	base = base.ReplaceAll(",","");
	base = base.ReplaceAll(".","");
	base.Remove(TString::kBoth,'_');

	int i=-1;
	TString uniqName("");

	do {
		++i;
		uniqName = base;
		uniqName+= "_";
		uniqName+= i;
	} while(usedNames.count(uniqName)>0);

	usedNames.insert(uniqName);
	
	return uniqName;
}

void RapidDecay::setupHistos() {
	for(unsigned int i=0; i<parts.size(); ++i) {
		if(nDaug[i] < 2) continue;
		// for each subdecay include:
		// parent mass
		// all 2- and 3-body daughter combinations
		// smeared versions of histograms
		TString baseName = names[i];
		TString histName;
		TString axisTitle;
		TH1F* hist(0);

		//for the mother mass
		double mmin = m[i];
		double mmax = m[i];
		int id = TMath::Abs(parts[i]);
		if(massdata.count(id)) {
			mmin = minmass[id];
			mmax = maxmass[id];
		}
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
		if(nDaug[i] > 2) {
			int first = firstDaug[i];

			for(int j=0; j<nDaug[i]; ++j) {
				for(int k=j+1; k<nDaug[i]; ++k) {
					mmin = m[first+j] + m[first+k];
					mmax = m[i];
					for(int daug=0; daug<nDaug[i]; ++daug) {
						if(daug==j || daug==k) continue;
						mmax -= m[first+daug];
					}
					axisTitle = "m(";
					axisTitle += j;
					axisTitle += k;
					axisTitle += ")";
					histName = baseName+"_M";
					histName += j;
					histName += k;
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
		if(nDaug[i] > 3) {
			int first = firstDaug[i];

			for(int j=0; j<nDaug[i]; ++j) {
				for(int k=j+1; k<nDaug[i]; ++k) {
					for(int l=k+1; l<nDaug[i]; ++l) {
						mmin = m[first+j] + m[first+k] + m[first+l];
						mmax = m[i];
						for(int daug=0; daug<nDaug[i]; ++daug) {
							if(daug==j || daug==k || daug==l) continue;
							mmax -= m[first+daug];
						}
						axisTitle = "m(";
						axisTitle += j;
						axisTitle += k;
						axisTitle += l;
						axisTitle += ")";
						histName = baseName+"_M";
						histName += j;
						histName += k;
						histName += l;
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
		if(nDaug[i] < 2) continue;
		// for each subdecay include:
		// parent mass
		// all 2- and 3-body daughter combinations
		// smeared versions of histograms

		//for the mother mass
		histos[iHist++]->Fill(p[i].M());
		histos[iHist++]->Fill(pSmeared[i].M());

		//2-body IMs
		if(nDaug[i] > 2) {
			int first = firstDaug[i];

			for(int j=0; j<nDaug[i]; ++j) {
				for(int k=j+1; k<nDaug[i]; ++k) {
					TLorentzVector pSum = p[first+j];
					pSum += p[first+k];
					histos[iHist++]->Fill(pSum.M());

					pSum = pSmeared[first+j];
					pSum += pSmeared[first+k];
					histos[iHist++]->Fill(pSum.M());
				}
			}
		}

		//3-body IMs
		if(nDaug[i] > 3) {
			int first = firstDaug[i];

			for(int j=0; j<nDaug[i]; ++j) {
				for(int k=j+1; k<nDaug[i]; ++k) {
					for(int l=k+1; l<nDaug[i]; ++l) {
						TLorentzVector pSum = p[first+j];
						pSum += p[first+k];
						pSum += p[first+l];
						histos[iHist++]->Fill(pSum.M());

						pSum = pSmeared[first+j];
						pSum += pSmeared[first+k];
						pSum += pSmeared[first+l];
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
		TString baseName = names[i];
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
		TLorentzVector mom = pSmeared[i];
		treeVars[varsPerPart*i+0] = parts[i];
		treeVars[varsPerPart*i+1] = mom.M();
		treeVars[varsPerPart*i+2] = mom.E();
		treeVars[varsPerPart*i+3] = mom.P();
		treeVars[varsPerPart*i+4] = mom.Px();
		treeVars[varsPerPart*i+5] = mom.Py();
		treeVars[varsPerPart*i+6] = mom.Pz();
		treeVars[varsPerPart*i+7] = mom.Pt();
		treeVars[varsPerPart*i+8] = mom.Eta();
		treeVars[varsPerPart*i+9] = mom.Phi();
		mom = p[i];
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
		int id = TMath::Abs(parts[i]);
		if(massdata.count(id)) {
			TString varName = "m";
			varName += id;
			m[i] = pick(massdata[id], rand, std::string(varName.Data()));
		}
	}
}

void RapidDecay::setupMasses() {//TODO store masses and widths of particles somewhere
	for(unsigned int i=0; i<parts.size(); ++i) {
		//check if it's already been loaded
		int id = TMath::Abs(parts[i]);
		if(massdata.count(id)) continue;

		switch(id) {
			case 213:
				setupRhoMass();
				break;
			case 323:
				setupKstMass();
				break;
			case 333:
				setupPhiMass();
				break;
			case 10441:
				setupChic0Mass();
				break;
			case 20443:
				setupChic1Mass();
				break;
			case 445:
				setupChic2Mass();
				break;
			default:
				break;
		}
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
		if(!genDecay()) continue;
		denomHisto->Fill(evalCustomParam(accRejParameter));
	}
	return denomHisto;
}

void RapidDecay::setupRhoMass() {
	RooRealVar m213("m213","m213",0.4, 1.5);
	RooGounarisSakurai* gs = createRhoPlus(m213);
	massdata[213] = gs->generate(RooArgSet(m213),100000);
	double mmin(0), mmax(0);
	massdata[213]->getRange(m213,mmin,mmax);
	minmass[213] = mmin;
	maxmass[213] = mmax;
}

void RapidDecay::setupKstMass() {
	RooRealVar m323("m323","m323",0.5, 1.5);
	RooRelBreitWigner* bw = createPhiMassPdf(m323);
	massdata[323] = bw->generate(RooArgSet(m323),100000);
	double mmin(0), mmax(0);
	massdata[323]->getRange(m323,mmin,mmax);
	minmass[323] = mmin;
	maxmass[323] = mmax;
}

void RapidDecay::setupPhiMass() {
	RooRealVar m333("m333","m333",0.6, 1.5);
	RooRelBreitWigner* bw = createPhiMassPdf(m333);
	massdata[333] = bw->generate(RooArgSet(m333),100000);
	double mmin(0), mmax(0);
	massdata[333]->getRange(m333,mmin,mmax);
	minmass[333] = mmin;
	maxmass[333] = mmax;
}

void RapidDecay::setupChic0Mass() {
	RooRealVar m10441("m10441","m10441",2.5, 5.0);
	RooRelBreitWigner* bw = createChi0MassPdf(m10441);
	massdata[10441] = bw->generate(RooArgSet(m10441),100000);
	double mmin(0), mmax(0);
	massdata[10441]->getRange(m10441,mmin,mmax);
	minmass[10441] = mmin;
	maxmass[10441] = mmax;
}

void RapidDecay::setupChic1Mass() {
	RooRealVar m20443("m20443","m20443",2.5, 5.0);
	RooRelBreitWigner* bw = createChi1MassPdf(m20443);
	massdata[20443] = bw->generate(RooArgSet(m20443),100000);
	double mmin(0), mmax(0);
	massdata[20443]->getRange(m20443,mmin,mmax);
	minmass[20443] = mmin;
	maxmass[20443] = mmax;
}

void RapidDecay::setupChic2Mass() {
	RooRealVar m445("m445","m445",2.5, 5.0);
	RooRelBreitWigner* bw = createChi2MassPdf(m445);
	massdata[445] = bw->generate(RooArgSet(m445),100000);
	double mmin(0), mmax(0);
	massdata[445]->getRange(m445,mmin,mmax);
	minmass[445] = mmin;
	maxmass[445] = mmax;
}

void RapidDecay::genParent() {
	double pt(0), eta(0), phi(rand.Uniform(0,2*TMath::Pi()));
	if(ptHisto)   pt = ptHisto->GetRandom();
	if(etaHisto) eta = etaHisto->GetRandom();
	p[0].SetPtEtaPhiM(pt,eta,phi,m[0]);
}

bool RapidDecay::genDecay() {
	int sumDaug(0);
	for(unsigned int i=0; i<parts.size(); ++i) {
		if(nDaug[i]>0) {
			TGenPhaseSpace event;
			if(!generateEvent(p[i], event, &m[1+sumDaug], nDaug[i], rand, maxgen)) {
				std::cout << "ERROR in RapidDecay::generate : generation failed." << std::endl;
				return false;
			}
			for(int j=0; j<nDaug[i]; ++j) {
				p[sumDaug+1+j] = *event.GetDecay(j);
			}

			//now increment the counters so we know where the next subdecay starts
			//++iDecay;
			sumDaug+=nDaug[i];
		}
	}

	return true;
}
