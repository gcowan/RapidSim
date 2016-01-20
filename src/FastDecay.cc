#include "FastDecay.h"

#include "TMath.h"

#include "RooRealVar.h"

#include <fstream>
#include <queue>

void FastDecay::loadParentKinematics(TH1F* pt, TH1F* eta) {
	std::cout << "INFO in FastDecay::loadParentKinematics : setting kinematics of the parent." << std::endl;
	ptHisto=pt;
	etaHisto=eta;
}

void FastDecay::setAcceptRejectHist(TH1F* hist, ParamType type, std::vector<int> particles) {
	if(accRejHisto) {
		std::cout << "WARNING in FastDecay::setAcceptRejectHist : accept/reject histogram already set. Original histogram will be used." << std::endl;
		return;
	}

	std::cout << "INFO in FastDecay::setAcceptRejectHist : setting the required kinematic distribution." << std::endl;
	accRejHisto = hist;
	CustomParameter param;
	param.type = type;
	param.particles = particles;
	accRejParameter = param;

	//correct the histogram to account for the the phasespace distribution
	accRejHisto->Divide(generateAccRejDenominator());
}

void FastDecay::addCustomParameter(TString name, ParamType type, std::vector<int> particles, bool truth, double min, double max) {
	std::cout << "INFO in FastDecay::addCustomParameter : adding custom parameter " << name << " = ";
	if(truth) std::cout << "TRUE";
	switch(type) {
		case FastDecay::M:
		       std::cout << "M";
		       break;
		case FastDecay::M2:
		       std::cout << "M2";
		       break;
		case FastDecay::MT:
		       std::cout << "MT";
		       break;
		case FastDecay::E:
			std::cout << "E";
		       break;
		case FastDecay::ET:
			std::cout << "ET";
		       break;
		case FastDecay::P:
		       std::cout << "P";
		       break;
		case FastDecay::PX:
		       std::cout << "PX";
		       break;
		case FastDecay::PY:
		       std::cout << "PY";
		       break;
		case FastDecay::PZ:
		       std::cout << "PZ";
		       break;
		case FastDecay::PT:
		       std::cout << "PT";
		       break;
		case FastDecay::ETA:
		       std::cout << "eta";
		       break;
		case FastDecay::PHI:
		       std::cout << "phi";
		       break;
		case FastDecay::RAPIDITY:
		       std::cout << "Rapidity";
		       break;
		case FastDecay::GAMMA:
		       std::cout << "gamma";
		       break;
		case FastDecay::BETA:
		       std::cout << "beta";
		       break;
	}
	std::cout << "(";
	for(int i=0; i<particles.size(); ++i) {
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

bool FastDecay::generate() {
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
		std::cout << "WARNING in FastDecay::generate : no events found with required kinematics." << std::endl;
		return false;
	}

	smearMomenta();
	fillHistos();
	if(tree) fillTree();

	return true;
}

void FastDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=p.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		if(nDaug[i] == 0) {
			//smear momentum of each detected particle
			pSmeared[i] = smearedVec(p[i],smearGraph,rand);
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
		
TLorentzVector FastDecay::getP(unsigned int i) {
	if(p.size() > i) return p[i];
	else {
		std::cout << "WARNING in FastDecay::getP : particle: " << i << "does not exist." << std::endl
			  << "                             Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}
		
TLorentzVector FastDecay::getPSmeared(unsigned int i) {
	if(pSmeared.size() > i) return pSmeared[i];
	else {
		std::cout << "WARNING in FastDecay::getPSmeared : particle: " << i << "does not exist." << std::endl
			  << "                                    Returning empty 4-vector..." << std::endl;
		return TLorentzVector();
	}
}

void FastDecay::saveHistos(TString fname) {
	std::cout << "INFO in FastDecay::saveHistos : saving histograms to file: " << fname << std::endl;
	TFile* histFile = new TFile(fname, "RECREATE");

	for(unsigned int i=0; i<histos.size(); ++i) {
		histos[i]->Write();
	}
	histFile->Close();
}

void FastDecay::saveTree(TString fname) {
	if(!tree) {
		std::cout << "INFO in FastDecay::saveTree : tree will be saved to file: " << fname << std::endl;
		std::cout << "                            : This will slow down generation." << std::endl;
		TFile* treeFile = new TFile(fname, "RECREATE");
		setupTree();
		tree->SetDirectory(treeFile);
	} else {
		std::cout << "WARNING in FastDecay::saveTree : tree has already been setup and will be saved to the original file." << std::endl;
	}
}

void FastDecay::loadDecay(TString filename) {
	std::cout << "INFO in FastDecay::loadDecay : loading decay descriptor from file: " << filename << std::endl;
	TString decayStr;
	std::queue<TString> decays;

	std::ifstream fin;
	fin.open(filename);
	decayStr.ReadLine(fin);
	fin.close();

	decays.push(decayStr);

	while(!decays.empty()) {

		decayStr = decays.front();
		std::cout << decayStr << std::endl;
		decays.pop();

		//first strip out any subdecays and add them to the queue
		while(decayStr.First('{')!=-1) {
			int start = decayStr.Index('{');
			int end = decayStr.Index('}',start);
			if(end < 0) {
				std::cout << "WARNING in FastDecay::loadDecay : malformed decay descriptor." << std::endl
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
			++count;

		}

		nDaug[theMother] = count;
		firstDaug[theMother] = parts.size() - count;

	}


	setupMasses();
	setupHistos();

	for(unsigned int i=0; i<parts.size(); ++i) {
		std::cout << i << "\t" << names[i] << "\t" << mother[i] << "\t" << parts[i] << "\t" << m[i] << "\t" << nDaug[i] << "\t" << firstDaug[i] << std::endl;
	}
}

TString FastDecay::getUniqName(TString base) {
	int i=-1;
	TString uniqName("");

	do {
		++i;
		uniqName = base;
		uniqName+= "_";
		uniqName+= i;
	} while(usedNames.count(uniqName)>0);

	//now sanitise
	uniqName = uniqName.ReplaceAll("+","p");
	uniqName = uniqName.ReplaceAll("-","m");
	uniqName = uniqName.ReplaceAll("*","st");
	uniqName = uniqName.ReplaceAll("(","_");
	uniqName = uniqName.ReplaceAll(")","_");
	uniqName = uniqName.ReplaceAll("[","_");
	uniqName = uniqName.ReplaceAll("]","_");
	uniqName = uniqName.ReplaceAll("<","_");
	uniqName = uniqName.ReplaceAll(">","_");
	uniqName = uniqName.ReplaceAll("{","_");
	uniqName = uniqName.ReplaceAll("}","_");
	uniqName = uniqName.ReplaceAll(" ","_");
	uniqName = uniqName.ReplaceAll("$","");
	uniqName = uniqName.ReplaceAll("%","");
	uniqName = uniqName.ReplaceAll("&","");
	uniqName = uniqName.ReplaceAll("/","");
	uniqName = uniqName.ReplaceAll(":","");
	uniqName = uniqName.ReplaceAll(";","");
	uniqName = uniqName.ReplaceAll("=","");
	uniqName = uniqName.ReplaceAll("\\","");
	uniqName = uniqName.ReplaceAll("^","");
	uniqName = uniqName.ReplaceAll("|","");
	uniqName = uniqName.ReplaceAll(",","");
	uniqName = uniqName.ReplaceAll(".","");
	uniqName.Remove(TString::kBoth,'_');

	usedNames.insert(uniqName);
	return uniqName;
}

void FastDecay::setupHistos() {
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

void FastDecay::fillHistos() {

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

void FastDecay::setupTree() {
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

void FastDecay::fillTree() {

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

void FastDecay::floatMasses() {
	for(unsigned int i=0; i<parts.size(); ++i) {
		int id = TMath::Abs(parts[i]);
		if(massdata.count(id)) {
			TString varName = "m";
			varName += id;
			m[i] = pick(massdata[id], rand, std::string(varName.Data()));
		}
	}
}

void FastDecay::setupMasses() {//TODO store masses and widths of particles somewhere
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

double FastDecay::evalCustomParam(int i) {
	CustomParameter param = customParams[i];
	evalCustomParam(param);
}

double FastDecay::evalCustomParam(CustomParameter param) {
	TLorentzVector mom;
	if(param.truth) {
		for(int i=0; i<param.particles.size(); ++i) {
			mom += getP(param.particles[i]);
		}
	} else {
		for(int i=0; i<param.particles.size(); ++i) {
			mom += getPSmeared(param.particles[i]);
		}
	}
	switch(param.type) {
		case FastDecay::M:
			return mom.M();
		case FastDecay::M2:
			return mom.M2();
		case FastDecay::MT:
			return mom.Mt();
		case FastDecay::E:
			return mom.E();
		case FastDecay::ET:
			return mom.Et();
		case FastDecay::P:
			return mom.P();
		case FastDecay::PX:
			return mom.Px();
		case FastDecay::PY:
			return mom.Py();
		case FastDecay::PZ:
			return mom.Pz();
		case FastDecay::PT:
			return mom.Pt();
		case FastDecay::ETA:
			return mom.Eta();
		case FastDecay::PHI:
			return mom.Phi();
		case FastDecay::RAPIDITY:
			return mom.Rapidity();
		case FastDecay::GAMMA:
			return mom.Gamma();
		case FastDecay::BETA:
			return mom.Beta();
	}

	return 0.;
}

bool FastDecay::runAcceptReject() {
	double val = evalCustomParam(accRejParameter);
	double score = accRejHisto->Interpolate(val);
	double max = accRejHisto->GetMaximum();
	if(score > rand.Uniform(max)) return true;
	return false;
}

TH1F* FastDecay::generateAccRejDenominator() {
	TH1F* denomHisto = dynamic_cast<TH1F*>(accRejHisto->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in FastDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay()) return false;
		denomHisto->Fill(evalCustomParam(accRejParameter));
	}
	return denomHisto;
}

void FastDecay::setupRhoMass() {
    RooRealVar m213("m213","m213",0.4, 1.5);
    RooGounarisSakurai* gs = createRhoPlus(m213);
    massdata[213] = gs->generate(RooArgSet(m213),100000);
    double mmin(0), mmax(0);
    massdata[213]->getRange(m213,mmin,mmax);
    minmass[213] = mmin;
    maxmass[213] = mmax;
}

void FastDecay::setupKstMass() {
    RooRealVar m323("m323","m323",0.5, 1.5);
    RooRelBreitWigner* bw = createPhiMassPdf(m323);
    massdata[323] = bw->generate(RooArgSet(m323),100000);
    double mmin(0), mmax(0);
    massdata[323]->getRange(m323,mmin,mmax);
    minmass[323] = mmin;
    maxmass[323] = mmax;
}

void FastDecay::setupPhiMass() {
    RooRealVar m333("m333","m333",0.6, 1.5);
    RooRelBreitWigner* bw = createPhiMassPdf(m333);
    massdata[333] = bw->generate(RooArgSet(m333),100000);
    double mmin(0), mmax(0);
    massdata[333]->getRange(m333,mmin,mmax);
    minmass[333] = mmin;
    maxmass[333] = mmax;
}

void FastDecay::setupChic0Mass() {
    RooRealVar m10441("m10441","m10441",2.5, 5.0);
    RooRelBreitWigner* bw = createChi0MassPdf(m10441);
    massdata[10441] = bw->generate(RooArgSet(m10441),100000);
    double mmin(0), mmax(0);
    massdata[10441]->getRange(m10441,mmin,mmax);
    minmass[10441] = mmin;
    maxmass[10441] = mmax;
}

void FastDecay::setupChic1Mass() {
    RooRealVar m20443("m20443","m20443",2.5, 5.0);
    RooRelBreitWigner* bw = createChi1MassPdf(m20443);
    massdata[20443] = bw->generate(RooArgSet(m20443),100000);
    double mmin(0), mmax(0);
    massdata[20443]->getRange(m20443,mmin,mmax);
    minmass[20443] = mmin;
    maxmass[20443] = mmax;
}

void FastDecay::setupChic2Mass() {
    RooRealVar m445("m445","m445",2.5, 5.0);
    RooRelBreitWigner* bw = createChi2MassPdf(m445);
    massdata[445] = bw->generate(RooArgSet(m445),100000);
    double mmin(0), mmax(0);
    massdata[445]->getRange(m445,mmin,mmax);
    minmass[445] = mmin;
    maxmass[445] = mmax;
}

void FastDecay::genParent() {
	double pt(0), eta(0), phi(rand.Uniform(0,2*TMath::Pi()));
	if(ptHisto)   pt = ptHisto->GetRandom();
	if(etaHisto) eta = etaHisto->GetRandom();
	p[0].SetPtEtaPhiM(pt,eta,phi,m[0]);
}

bool FastDecay::genDecay() {
	int sumDaug(0);
	for(unsigned int i=0; i<parts.size(); ++i) {
		if(nDaug[i]>0) {
			TGenPhaseSpace event;
			if(!generateEvent(p[i], event, &m[1+sumDaug], nDaug[i], rand, maxgen)) {
				std::cout << "ERROR in FastDecay::generate : generation failed." << std::endl;
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
