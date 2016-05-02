#include "RapidDecay.h"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParam.h"
#include "RapidParticle.h"
#include "RapidParticleData.h"

void RapidDecay::setParentKinematics(TH1* ptHisto, TH1* etaHisto) {
	std::cout << "INFO in RapidDecay::setParentKinematics : setting kinematics of the parent." << std::endl;
	ptHisto_=ptHisto;
	etaHisto_=etaHisto;
}

void RapidDecay::setAcceptRejectHist(TH1* histo, RapidParam* param) {
	accRejParameterX_ = param;
	accRejParameterY_ = 0;
	accRejHisto_      = histo;

	//correct the histogram to account for the the phasespace distribution
	TH1* denom = generateAccRejDenominator1D();
	accRejHisto_->Divide(denom);
	delete denom;
}

void RapidDecay::setAcceptRejectHist(TH1* histo, RapidParam* paramX, RapidParam* paramY) {
	accRejParameterX_ = paramX;
	accRejParameterY_ = paramY;
	accRejHisto_      = histo;

	//correct the histogram to account for the the phasespace distribution
	TH2* denom = generateAccRejDenominator2D();
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

	return true;
}

void RapidDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=parts_.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		parts_[i]->smearMomentum();
	}

}

void RapidDecay::setup() {
	setupMasses();

	std::cout << "INFO in RapidDecay::setup : Particle summary follows:" << std::endl;
	printf("index\tlabel\t\t   ID\t\tmass (GeV/c^2)\tmother\t\t# daughters\tdaughters\n");
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->print(i);
	}
}

void RapidDecay::floatMasses() {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->floatMass();
	}
}

void RapidDecay::setupMasses() {
	RapidParticleData* particleData = RapidParticleData::getInstance();

	for(unsigned int i=0; i<parts_.size(); ++i) {
		particleData->setupMass(parts_[i]);
	}
}

bool RapidDecay::runAcceptReject() {
	if(accRejParameterY_) return runAcceptReject2D();
	else return runAcceptReject1D();
}

bool RapidDecay::runAcceptReject1D() {
	double val = accRejParameterX_->eval();
	int bin = accRejHisto_->FindBin(val);

	double score(0.);
	if(!accRejHisto_->IsBinOverflow(bin) && !accRejHisto_->IsBinUnderflow(bin)) {
		score = accRejHisto_->Interpolate(val);
	}
	double max = accRejHisto_->GetMaximum();
	if(score > gRandom->Uniform(max)) return true;
	return false;
}

bool RapidDecay::runAcceptReject2D() {
	double valX = accRejParameterX_->eval();
	double valY = accRejParameterY_->eval();
	int bin = accRejHisto_->FindBin(valX,valY);

	double score(0.);
	if(!accRejHisto_->IsBinOverflow(bin) && !accRejHisto_->IsBinUnderflow(bin)) {
		score = accRejHisto_->Interpolate(valX,valY);
	}
	double max = accRejHisto_->GetMaximum();
	if(score > gRandom->Uniform(max)) return true;
	return false;
}

TH1* RapidDecay::generateAccRejDenominator1D() {
	TH1* denomHisto = dynamic_cast<TH1*>(accRejHisto_->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in RapidDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay(true)) continue;
		denomHisto->Fill(accRejParameterX_->eval());
	}
	return denomHisto;
}

TH2* RapidDecay::generateAccRejDenominator2D() {
	TH2* denomHisto = dynamic_cast<TH2*>(accRejHisto_->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in RapidDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay(true)) continue;
		denomHisto->Fill(accRejParameterX_->eval(), accRejParameterY_->eval());
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
			// check decay kinematics valid
			if(!decay_.SetDecay(part->getP(), part->nDaughters(), part->daughterMasses())) {
				std::cout << "ERROR in RapidDecay::genDecay : decay of " << part->name() << " is kinematically forbidden." << std::endl;
				return false;
			}

			// make an event
			if(acceptAny) {
				decay_.Generate();
			} else {
				int nGen(0);
				bool accept(false);
				while (nGen < maxgen_ && accept == false){
					accept = decay_.Generate() > gRandom->Uniform();
					++nGen;
				} // while

				if(!accept) {
					std::cout << "ERROR in RapidDecay::genDecay : rejected all " << maxgen_ << " attempts to decay " << part->name() << "." << std::endl;
					return false;
				}
			}

			int j=0;
			for(RapidParticle* jDaug=part->daughter(0); jDaug!=0; jDaug=jDaug->next()) {
				jDaug->setP(*decay_.GetDecay(j++));
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
