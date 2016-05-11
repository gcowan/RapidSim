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

void RapidHistWriter::setup(bool saveTree) {
	setupHistos();
	if(saveTree) setupTree(); //called after setupHistos so we know the number of parameters
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
	std::vector<RapidParam*>::iterator itParam;

	TString baseName;
	TString histName;
	TString axisTitle;
	TH1F* hist(0);

	double min(0), max(0);

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		baseName = part->name();

		if(part->nDaughters() < 2) {
			for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
				RapidParam* param = *itParam;
				histName = baseName+"_"+param->name();
				axisTitle = param->name() + "(" + baseName + ")";
				param->getMinMax(part,min,max);
				hist = new TH1F(histName, "", 100, min, max);
				hist->GetXaxis()->SetTitle(axisTitle);
				histos_.push_back(hist);
			}
		} else {
			for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
				RapidParam* param = *itParam;
				histName = baseName+"_"+param->name();
				axisTitle = param->name() + "(" + baseName + ")";
				param->getMinMax(part,min,max);
				hist = new TH1F(histName, "", 100, min, max);
				hist->GetXaxis()->SetTitle(axisTitle);
				histos_.push_back(hist);
			}

			//2-body IMs
			if(part->nDaughters() > 2) {
				RapidParticle* jDaug = part->daughter(0);

				for( ; jDaug!=0; jDaug=jDaug->next()) {
					for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
						for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
							RapidParam* param = *itParam;
							histName  = jDaug->name()+"_";
							histName += kDaug->name()+"_";
							histName += param->name();
							axisTitle = param->name() + "(" + jDaug->name() + "_" + kDaug->name() + ")";
							param->getMinMax(jDaug,kDaug,min,max);
							hist = new TH1F(histName, "", 100, min, max);
							hist->GetXaxis()->SetTitle(axisTitle);
							histos_.push_back(hist);
						}
					}
				}
			}

			//3-body IMs
			if(part->nDaughters() > 3) {
				RapidParticle* jDaug = part->daughter(0);

				for( ; jDaug!=0; jDaug=jDaug->next()) {
					for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
						for(RapidParticle* lDaug=kDaug->next(); lDaug!=0; lDaug=lDaug->next()) {
							for(itParam=paramsThreeBody_.begin(); itParam!=paramsThreeBody_.end(); ++itParam) {
								RapidParam* param = *itParam;
								histName  = jDaug->name()+"_";
								histName += kDaug->name()+"_";
								histName += lDaug->name()+"_";
								histName += param->name();
								axisTitle = param->name() + "(" + jDaug->name() + "_" + kDaug->name() + "_" + lDaug->name() + ")";
								param->getMinMax(jDaug,kDaug,lDaug,min,max);
								hist = new TH1F(histName, "", 100, min, max);
								hist->GetXaxis()->SetTitle(axisTitle);
								histos_.push_back(hist);
							}
						}

					}
				}
			}
		}
	}

	for( itParam=params_.begin(); itParam!= params_.end(); ++itParam) {
		RapidParam* param = (*itParam);
		TH1F* hist = new TH1F(param->name(), "", 100, param->min(), param->max());
		histos_.push_back(hist);
	}

}

void RapidHistWriter::setupTree() {
	std::cout << "INFO in RapidHistWriter::setupTree : tree will be saved to file: " << name_ << "_tree.root" << std::endl;
	std::cout << "                                     This will slow down generation." << std::endl;

	treeFile_ = new TFile(name_+"_tree.root", "RECREATE");
	tree_ = new TTree("DecayTree","DecayTree");
	tree_->SetDirectory(treeFile_);

	treeVars_ = std::vector<double>(histos_.size(), 0);
	for(unsigned int i=0; i<histos_.size(); ++i) {
		tree_->Branch(histos_[i]->GetName(), &treeVars_[i]);
	}

}

void RapidHistWriter::fillHistos() {

	int iHist(0);
	std::vector<RapidParam*>::iterator itParam;
	TLorentzVector pSumTruth;
	TLorentzVector pSum;

	double val;

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];

		if(part->nDaughters() < 2) {
			for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
				RapidParam* param = *itParam;
				if(param->truth()) {
					val = param->eval(part->getP());
				} else {
					val = param->eval(part->getPSmeared());
				}
				histos_[iHist++]->Fill(val);
			}
		} else {
			for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
				RapidParam* param = *itParam;
				if(param->truth()) {
					val = param->eval(part->getP());
				} else {
					val = param->eval(part->getPSmeared());
				}
				histos_[iHist++]->Fill(val);
			}

			//2-body IMs
			if(part->nDaughters() > 2) {
				RapidParticle* jDaug = part->daughter(0);

				for( ; jDaug!=0; jDaug=jDaug->next()) {
					for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
						pSumTruth = jDaug->getP() + kDaug->getP();
						pSum = jDaug->getPSmeared() + kDaug->getPSmeared();

						for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
							RapidParam* param = *itParam;
							if(param->truth()) {
								val = param->eval(pSumTruth);
							} else {
								val = param->eval(pSum);
							}
							histos_[iHist++]->Fill(val);
						}
					}
				}
			}

			//3-body IMs
			if(part->nDaughters() > 3) {
				RapidParticle* jDaug = part->daughter(0);

				for( ; jDaug!=0; jDaug=jDaug->next()) {
					for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
						for(RapidParticle* lDaug=kDaug->next(); lDaug!=0; lDaug=lDaug->next()) {
							pSumTruth = jDaug->getP() + kDaug->getP() + lDaug->getP();
							pSum = jDaug->getPSmeared() + kDaug->getPSmeared() + lDaug->getPSmeared();

							for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
								RapidParam* param = *itParam;
								if(param->truth()) {
									val = param->eval(pSumTruth);
								} else {
									val = param->eval(pSum);
								}
								histos_[iHist++]->Fill(val);
							}
						}
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
	int iVar(0);
	std::vector<RapidParam*>::iterator itParam;
	TLorentzVector pSumTruth;
	TLorentzVector pSum;

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];

		if(part->nDaughters() < 2) {
			for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
				RapidParam* param = *itParam;
				if(param->truth()) {
					treeVars_[iVar++] = param->eval(part->getP());
				} else {
					treeVars_[iVar++] = param->eval(part->getPSmeared());
				}
			}
		} else {
			for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
				RapidParam* param = *itParam;
				if(param->truth()) {
					treeVars_[iVar++] = param->eval(part->getP());
				} else {
					treeVars_[iVar++] = param->eval(part->getPSmeared());
				}
			}

			//2-body IMs
			if(part->nDaughters() > 2) {
				RapidParticle* jDaug = part->daughter(0);

				for( ; jDaug!=0; jDaug=jDaug->next()) {
					for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
						pSumTruth = jDaug->getP() + kDaug->getP();
						pSum = jDaug->getPSmeared() + kDaug->getPSmeared();

						for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
							RapidParam* param = *itParam;
							if(param->truth()) {
								treeVars_[iVar++] = param->eval(pSumTruth);
							} else {
								treeVars_[iVar++] = param->eval(pSum);
							}
						}
					}
				}
			}

			//3-body IMs
			if(part->nDaughters() > 3) {
				RapidParticle* jDaug = part->daughter(0);

				for( ; jDaug!=0; jDaug=jDaug->next()) {
					for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
						for(RapidParticle* lDaug=kDaug->next(); lDaug!=0; lDaug=lDaug->next()) {
							pSumTruth = jDaug->getP() + kDaug->getP() + lDaug->getP();
							pSum = jDaug->getPSmeared() + kDaug->getPSmeared() + lDaug->getPSmeared();

							for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
								RapidParam* param = *itParam;
								if(param->truth()) {
									treeVars_[iVar++] = param->eval(pSumTruth);
								} else {
									treeVars_[iVar++] = param->eval(pSum);
								}
							}
						}
					}
				}
			}
		}
	}

	std::vector<RapidParam*>::iterator it = params_.begin();
	for( ; it!= params_.end(); ++it) {
		treeVars_[iVar++] = (*it)->eval();
	}

	tree_->Fill();
}
