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

	if(tree_) tree_->Fill();
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

	TString histName;
	TH1F* hist(0);
	double min(0), max(0);

	// Following code could be made better by passing in the list of params
	for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
		histName = (*itParam)->name() + suffix;
		(*itParam)->getMinMax(min,max,true);
		hist = new TH1F(histName, "", 100, min, max);
		histos_.push_back(hist);
	}

	for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
		histName = (*itParam)->name() + suffix;
		(*itParam)->getMinMax(min,max,true);
		hist = new TH1F(histName, "", 100, min, max);
		histos_.push_back(hist);
	}

	for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
		histName = (*itParam)->name() + suffix;
		(*itParam)->getMinMax(min,max,true);
		hist = new TH1F(histName, "", 100, min, max);
		histos_.push_back(hist);
	}

	for(itParam=paramsThreeBody_.begin(); itParam!=paramsThreeBody_.end(); ++itParam) {
		histName = (*itParam)->name() + suffix;
		(*itParam)->getMinMax(min,max,true);
		hist = new TH1F(histName, "", 100, min, max);
		histos_.push_back(hist);
	}

	for( itParam=params_.begin(); itParam!= params_.end(); ++itParam) {
		histName = (*itParam)->name() + suffix;
		(*itParam)->getMinMax(min,max,true);
		TH1F* hist = new TH1F(histName, "", 100, min, max);
		histos_.push_back(hist);
	}
}

void RapidHistWriter::setupTree() {
	std::cout << "INFO in RapidHistWriter::setupTree : tree will be saved to file: " << name_ << "_tree.root" << std::endl;
	std::cout << "                                     This will slow down generation." << std::endl;

	treeFile_ = new TFile(name_+"_tree.root", "RECREATE");
	tree_ = new TTree("DecayTree","DecayTree");
	tree_->SetDirectory(treeFile_);

    tree_->Branch("nEvent",&nevent_);
	for(unsigned int i=0; i<histos_.size(); ++i) {
		tree_->Branch(histos_[i]->GetName(), &vars_[i]);
	}

}

unsigned int RapidHistWriter::fillSingleHypothesis(unsigned int offset) {

	std::vector<RapidParam*>::iterator itParam;

	// Following code could be made a bit nicer by passing the params list
	// Note that the ordering here must be the same as in setupSingleHypothesis
	for(itParam=paramsDecaying_.begin(); itParam!=paramsDecaying_.end(); ++itParam) {
		RapidParam* param = *itParam;
		vars_[offset] = param->eval();
		histos_[offset]->Fill(vars_[offset]);
		++offset;
	}
	for(itParam=paramsStable_.begin(); itParam!=paramsStable_.end(); ++itParam) {
		RapidParam* param = *itParam;
		vars_[offset] = param->eval();
		histos_[offset]->Fill(vars_[offset]);
		++offset;
	}
	//2-body IMs
	for(itParam=paramsTwoBody_.begin(); itParam!=paramsTwoBody_.end(); ++itParam) {
		RapidParam* param = *itParam;
		vars_[offset] = param->eval();
		histos_[offset]->Fill(vars_[offset]);
		++offset;
	}

	//3-body IMs
	for(itParam=paramsThreeBody_.begin(); itParam!=paramsThreeBody_.end(); ++itParam) {
		RapidParam* param = *itParam;
		vars_[offset] = param->eval();
		histos_[offset]->Fill(vars_[offset]);
		++offset;
	}

	std::vector<RapidParam*>::iterator it = params_.begin();
	for( ; it!= params_.end(); ++it) {
		vars_[offset] = (*it)->eval();
		histos_[offset]->Fill(vars_[offset]);
		++offset;
	}

	return offset;
}
