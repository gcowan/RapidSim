#include "RapidAcceptance.h"

#include <iostream>

#include "TMath.h"

#include "RapidCut.h"
#include "RapidParticle.h"

RapidAcceptance::AcceptanceType RapidAcceptance::typeFromString(TString str) {
	if(str=="Any") {
		return RapidAcceptance::ANY;
	} else if(str=="MotherIn") {
		return RapidAcceptance::MOTHERIN;
	} else if(str=="AllIn") {
		return RapidAcceptance::ALLIN;
	} else if(str=="AllDownstream") {
		return RapidAcceptance::ALLDOWNSTREAM;
	} else {
		std::cout << "WARNING in RapidAcceptance::typeFromString : unknown type name " << str << "." << std::endl
			  << "                                             returning \"accept any event\" type." << std::endl;
		return RapidAcceptance::ANY;
	}
}

bool RapidAcceptance::isSelected() {
	if(!inAcceptance()) {
		return false;
	}

	std::vector<RapidCut*>::iterator it = cuts_.begin();
	for( ; it!= cuts_.end(); ++it) {
		if(!(*it)->passCut()) {
			return false;
		}
	}

	return true;
}

bool RapidAcceptance::inAcceptance() {
	switch(type_) {
		case MOTHERIN:
			return motherInAcceptance();
		case ALLIN:
			return allInAcceptance();
		case ALLDOWNSTREAM:
			return allInDownstream();
		case ANY:
		default:
			return true;
	}
}

void RapidAcceptance::setup(std::vector<RapidParticle*> parts) {

	//keep top particle
	if(!parts.empty()) top_ = parts[0];

	std::vector<RapidParticle*>::iterator it = parts.begin();
	for( ; it!= parts.end(); ++it) {
		RapidParticle* part = *it;
		//keep all stable particles
		if(part->stable()) parts_.push_back(part);

	}
}

bool RapidAcceptance::motherInAcceptance() {
	return partInAcceptance(top_);
}

bool RapidAcceptance::allInAcceptance() {
	std::vector<RapidParticle*>::iterator it = parts_.begin();
	for( ; it!= parts_.end(); ++it) {
		RapidParticle* part = *it;
		if(!partInAcceptance(part)) return false;
	}
	return true;
}

bool RapidAcceptance::allInDownstream() {
	if(!allInAcceptance()) return false;

	std::vector<RapidParticle*>::iterator it = parts_.begin();
	for( ; it!= parts_.end(); ++it) {
		RapidParticle* part = *it;
		if(!partInDownstream(part)) return false;
	}
	return true;
}

bool RapidAcceptance::partInAcceptance(RapidParticle* part) {

	TLorentzVector vec = part->getP();

	if (TMath::Abs(vec.Px()/vec.Pz()) > 0.3) return false;
	if (TMath::Abs(vec.Py()/vec.Pz()) > 0.25) return false;
	if (sqrt(pow(vec.Px()/vec.Pz(),2) + pow(vec.Py()/vec.Pz(),2)) <0.01) return false;

	return true;
}

bool RapidAcceptance::partInDownstream(RapidParticle* part) {

	TLorentzVector vec = part->getP();
	double charge = part->charge();
	TLorentzVector newvec = magnetKick(vec,charge);

	// position at magnet centre
	double xMag = zC_*vec.Px()/vec.Pz();
	double yMag = zC_*vec.Py()/vec.Pz();

	double xTracker = xMag + (newvec.Px()*(zTracker_ - zC_)/newvec.Pz());
	double yTracker = yMag + (newvec.Py()*(zTracker_ - zC_)/newvec.Pz());

	bool inX = TMath::Abs(xTracker) < xSizeTracker_ && TMath::Abs(xTracker) > xMinTracker_;
	bool inY = TMath::Abs(yTracker) < ySizeTracker_ && TMath::Abs(yTracker) > yMinTracker_;

	return inX && inY;
}

TLorentzVector RapidAcceptance::magnetKick(TLorentzVector& vec, double charge){

	// kick
	double p = vec.P();
	double px = vec.Px() + ptkick_*charge;
	double py = vec.Py();
	double pz = sqrt(p*p - px*px - py*py);

	TLorentzVector kicked;
	kicked.SetXYZM(px,py,pz,vec.M());

	return kicked;
}
