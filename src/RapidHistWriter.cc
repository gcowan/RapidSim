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
	//get the subset of particles that have more than one mass hypothesis
	std::vector<RapidParticle*>::const_iterator it = parts_.begin();
	for( ; it!=parts_.end(); ++it) {
		RapidParticle* part = *it;
		if(part->nMassHypotheses()>1) {
			altHypothesisParts_.push_back(part);
		}
	}

	setupHistos();
	if(saveTree) setupTree(); //called after setupHistos so we know the number of parameters
}

void RapidHistWriter::fill() {
	unsigned int offset(0);

	//first fill for the default hypothesis
	offset = fillSingleHypothesis(offset);

	RapidParticle *part1(0), *part2(0);

	//now loop over all single mis-IDs
	std::vector<RapidParticle*>::const_iterator it1 = altHypothesisParts_.begin();
	for( ; it1!=altHypothesisParts_.end(); ++it1) {
		part1 = *it1;
		for(unsigned int i=1; i<part1->nMassHypotheses(); ++i) {
			part1->setMassHypothesis(i);
			offset = fillSingleHypothesis(offset);
		}
		part1->setMassHypothesis(0);
	}

	//now loop over all double mis-IDs
	it1 = altHypothesisParts_.begin();
	for( ; it1!=altHypothesisParts_.end(); ++it1) {
		part1 = *it1;
		std::vector<RapidParticle*>::const_iterator it2 = it1+1;
		for( ; it2!=altHypothesisParts_.end(); ++it2) {
			part2 = *it2;
			for(unsigned int i=1; i<part1->nMassHypotheses(); ++i) {
				part1->setMassHypothesis(i);
				for(unsigned int j=1; j<part2->nMassHypotheses(); ++j) {
					part2->setMassHypothesis(j);
					offset = fillSingleHypothesis(offset);
				}
				part2->setMassHypothesis(0);
			}
			part1->setMassHypothesis(0);
		}
	}
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
	//first setup histograms for default mass hypothesis
	setupSingleHypothesis();

	RapidParticle *part1(0), *part2(0);

	//now loop over all single mis-IDs
	std::vector<RapidParticle*>::const_iterator it1 = altHypothesisParts_.begin();
	for( ; it1!=altHypothesisParts_.end(); ++it1) {
		part1 = *it1;
		for(unsigned int i=1; i<part1->nMassHypotheses(); ++i) {
			part1->setMassHypothesis(i);
			TString suffix("_"); suffix+=part1->name() + "2" + part1->massHypothesisName();
			setupSingleHypothesis(suffix);
		}
		part1->setMassHypothesis(0);
	}

	//now loop over all double mis-IDs
	it1 = altHypothesisParts_.begin();
	for( ; it1!=altHypothesisParts_.end(); ++it1) {
		part1 = *it1;
		std::vector<RapidParticle*>::const_iterator it2 = it1+1;
		for( ; it2!=altHypothesisParts_.end(); ++it2) {
			part2 = *it2;
			for(unsigned int i=1; i<part1->nMassHypotheses(); ++i) {
				part1->setMassHypothesis(i);
				for(unsigned int j=1; j<part2->nMassHypotheses(); ++j) {
					part2->setMassHypothesis(j);
					TString suffix("_"); suffix+=part1->name() + "2" + part1->massHypothesisName();
					suffix+="_"; suffix+=part2->name() + "2" + part2->massHypothesisName();
					setupSingleHypothesis(suffix);
				}
				part2->setMassHypothesis(0);
			}
			part1->setMassHypothesis(0);
		}
	}

	vars_ = std::vector<double>(histos_.size(), 0);
}


void RapidHistWriter::setupSingleHypothesis(TString suffix) {
	std::vector<RapidParam*>::iterator itParam;

	TString baseName;
	TString histName;
	TString axisTitle;
	TH1F* hist(0);

	double min(0), max(0);

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		baseName = part->name();

		if(part->nDaughters()==0) {
			for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
				RapidParam* param = *itParam;
				histName  = baseName+"_"+param->name();
				histName += suffix;
				axisTitle = param->name() + "(" + baseName + ")";
				param->getMinMax(part,min,max);
				hist = new TH1F(histName, "", 100, min, max);
				hist->GetXaxis()->SetTitle(axisTitle);
				histos_.push_back(hist);
			}
		} else {
			for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
				RapidParam* param = *itParam;
				histName  = baseName+"_"+param->name();
				histName += suffix;
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
							histName += suffix;
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
								histName += suffix;
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
		param->getMinMax(min,max,true);
		TH1F* hist = new TH1F(param->name()+suffix, "", 100, min, max);
		histos_.push_back(hist);
	}
}

void RapidHistWriter::setupTree() {
	std::cout << "INFO in RapidHistWriter::setupTree : tree will be saved to file: " << name_ << "_tree.root" << std::endl;
	std::cout << "                                     This will slow down generation." << std::endl;

	treeFile_ = new TFile(name_+"_tree.root", "RECREATE");
	tree_ = new TTree("DecayTree","DecayTree");
	tree_->SetDirectory(treeFile_);

	for(unsigned int i=0; i<histos_.size(); ++i) {
		tree_->Branch(histos_[i]->GetName(), &vars_[i]);
	}

}

unsigned int RapidHistWriter::fillSingleHypothesis(unsigned int offset) {

	std::vector<RapidParam*>::iterator itParam;
	TLorentzVector pSumTruth;
	TLorentzVector pSum;

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];

		if(part->nDaughters() < 2) {
			for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
				RapidParam* param = *itParam;
				if(param->truth()) {
					vars_[offset] = param->eval(part->getP());
				} else {
					vars_[offset] = param->eval(part->getPSmeared());
				}

				histos_[offset]->Fill(vars_[offset]);
				++offset;
			}
		} else {
			for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
				RapidParam* param = *itParam;
				if(param->truth()) {
					vars_[offset] = param->eval(part->getP());
				} else {
					vars_[offset] = param->eval(part->getPSmeared());
				}
				histos_[offset]->Fill(vars_[offset]);
				++offset;
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
								vars_[offset] = param->eval(pSumTruth);
							} else {
								vars_[offset] = param->eval(pSum);
							}
							histos_[offset]->Fill(vars_[offset]);
							++offset;
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
									vars_[offset] = param->eval(pSumTruth);
								} else {
									vars_[offset] = param->eval(pSum);
								}
								histos_[offset]->Fill(vars_[offset]);
								++offset;
							}
						}
					}
				}
			}
		}
	}

	std::vector<RapidParam*>::iterator it = params_.begin();
	for( ; it!= params_.end(); ++it) {
		vars_[offset] = (*it)->eval();
		histos_[offset]->Fill(vars_[offset]);
		++offset;
	}

	if(tree_) tree_->Fill();

	return offset;
}
