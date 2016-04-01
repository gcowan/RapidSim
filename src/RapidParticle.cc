#include "RapidParticle.h"

#include "TMath.h"

#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"

void RapidParticle::addDaughter(RapidParticle* part) {
	if(!daughters_.empty()) {
		daughters_[daughters_.size()-1]->next_ = part;
	}
	daughters_.push_back(part);
	daughterMasses_.push_back(part->mass());
}

void RapidParticle::smearMomentum() {
	if(nDaughters() == 0) {
		//smear momentum
		if(invisible_) {
			pSmeared_.SetXYZM(0.,0.,0.,0.);
		} else if(momSmear_) {
			pSmeared_ = momSmear_->smearMomentum(p_);
		} else {
			pSmeared_ = p_;
		}
	} else {
		//reconstruct mothers from their daughters
		pSmeared_ = TLorentzVector();
		RapidParticle* daug = daughter(0);
		for( ; daug!=0; daug = daug->next()) {
			pSmeared_ += daug->pSmeared_;
		}
	}
}

RapidParticle* RapidParticle::daughter(unsigned int i) {
	if(daughters_.size()>i+1) {
		return daughters_[i];
	}
	return 0;
}

void RapidParticle::print(int index) 
{
	TString mname = "---";
	if(mother_) mname = mother_->name_;
	TString dname = "";
	RapidParticle* daug = daughter(0);
	if(!daug) dname = "---";
	for( ; daug!=0; daug=daug->next()) {
		dname += daug->name_;
		if(daug->next()) dname+=", ";
	}

	printf("%3d\t%-15s\t%6d\t\t%.6f\t%-15s\t%2d\t\t%-15s\n", index, name_.Data(), id_, mass_, mname.Data(), nDaughters(), dname.Data());
}
