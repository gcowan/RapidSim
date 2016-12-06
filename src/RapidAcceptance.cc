#include "RapidAcceptance.h"

#include <iostream>

#include "TMath.h"

#include "RapidCut.h"
#include "RapidParticle.h"

RapidAcceptance::AcceptanceType RapidAcceptance::typeFromString(TString str) {
	if(str=="Any") {
		return RapidAcceptance::ANY;
	} else if(str=="ParentIn") {
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

RapidAcceptance::DetectorType RapidAcceptance::detectorFromString(TString str) {
	if(str=="LHCb") {
		return RapidAcceptance::LHCB;
	} else {
		std::cout << "WARNING in RapidAcceptance::detectorFromString : unknown detector name " << str << "." << std::endl
			  << "                                                 returning \"LHCb\"." << std::endl;
		return RapidAcceptance::LHCB;
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
