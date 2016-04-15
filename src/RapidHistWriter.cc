#include "RapidHistWriter.h"

#include "RapidParam.h"
#include "RapidParticle.h"

RapidHistWriter::~RapidHistWriter() {
	if(tree_) {
		tree_->AutoSave();
		treeFile_->Close();
	}
	while(!histos_.empty()) {
		delete histos_[histos_.size()-1];
		histos_.pop_back();
	}
}

void RapidHistWriter::setup(std::vector<RapidParticle*> parts, std::vector<RapidParam*> params, TString name, bool saveTree) {
	parts_ = parts;
	params_ = params;
	name_ = name;

	setupHistos();
	if(saveTree) setupTree();
}

void RapidHistWriter::fill() {
	fillHistos();
	if(tree_) fillTree();
}

void RapidHistWriter::save() {
	std::cout << "INFO in RapidHistWriter::save : saving histograms to file: " << name_ << "_hists.root" << std::endl;
	TFile* histFile = new TFile(name_+"_hists.root", "RECREATE");

	for(unsigned int i=0; i<histos_.size(); ++i) {
		histos_[i]->Write();
	}
	histFile->Close();

	if(tree_) {
		tree_->AutoSave();
	}
}

void RapidHistWriter::setupHistos() {
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

	std::vector<RapidParam*>::iterator it = params_.begin();
	for( ; it!= params_.end(); ++it) {
		RapidParam* param = (*it);
		TH1F* hist = new TH1F(param->name(), "", 100, param->min(), param->max());
		histos_.push_back(hist);
	}

}

void RapidHistWriter::setupTree() {
	std::cout << "INFO in RapidHistWriter::setupTree : tree will be saved to file: " << name_ << "_tree.root" << std::endl;
	std::cout << "                                   : This will slow down generation." << std::endl;

	varsPerPart_ = 19;
	treeFile_ = new TFile(name_+"_tree.root", "RECREATE");
	tree_ = new TTree("DecayTree","DecayTree");
	tree_->SetDirectory(treeFile_);
	treeVars_ = std::vector<double>(parts_.size()*varsPerPart_ + params_.size(), 0);

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
	for(unsigned int i=0; i<params_.size(); ++i) {
		tree_->Branch(params_[i]->name(), &treeVars_[varsPerPart_*parts_.size() + i]);
	}
}

void RapidHistWriter::fillHistos() {

	int iHist(0);

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		if(part->nDaughters() < 2) continue;
		// for each subdecay include:
		// parent mass
		// all 2- and 3-body daughter combinations
		// smeared versions of histograms

		//for the mother mass
		histos_[iHist++]->Fill(parts_[i]->getP().M());
		histos_[iHist++]->Fill(parts_[i]->getPSmeared().M());

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

	std::vector<RapidParam*>::iterator it = params_.begin();
	for( ; it!= params_.end(); ++it) {
		histos_[iHist++]->Fill((*it)->eval());
	}
}

void RapidHistWriter::fillTree() {

	for(unsigned int i=0; i<parts_.size(); ++i) {
		TLorentzVector mom = parts_[i]->getPSmeared();
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
		mom = parts_[i]->getP();
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

	for(unsigned int i=0; i<params_.size(); ++i) {
		treeVars_[varsPerPart_*parts_.size() + i] = params_[i]->eval();
	}

	tree_->Fill();
}
