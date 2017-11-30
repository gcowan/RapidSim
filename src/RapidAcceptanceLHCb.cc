#include "RapidAcceptanceLHCb.h"

#include <iostream>

#include "TMath.h"

#include "RapidParticle.h"

void RapidAcceptanceLHCb::getDefaultPtRange(double& min, double& max) {
	std::cout << "INFO in RapidAcceptanceLHCb::getDefaultPtRange : Getting pT range for LHCb geometry." << std::endl;
	std::cout << "                                                 Range is 0 - 100 GeV." << std::endl;
	min=0.;
	max=100.;
}
void RapidAcceptanceLHCb::getDefaultEtaRange(double& min, double& max) {
	std::cout << "INFO in RapidAcceptanceLHCb::getDefaultEtaRange : Getting eta range for LHCb geometry." << std::endl;
	std::cout << "                                                  Range is 1.0 - 6.0." << std::endl;
	min=1.;
	max=6.;
}

bool RapidAcceptanceLHCb::partInAcceptance(RapidParticle* part) {

	if(part->invisible()) return true;

	TLorentzVector vec = part->getP();

	if (TMath::Abs(vec.Px()/vec.Pz()) > 0.3) return false;
	if (TMath::Abs(vec.Py()/vec.Pz()) > 0.25) return false;
	if (sqrt(pow(vec.Px()/vec.Pz(),2) + pow(vec.Py()/vec.Pz(),2)) <0.01) return false;

	return true;
}

bool RapidAcceptanceLHCb::partInDownstream(RapidParticle* part) {

	if(part->invisible()) return true;

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

TLorentzVector RapidAcceptanceLHCb::magnetKick(TLorentzVector& vec, double charge){

	// kick
	double p = vec.P();
	double px = vec.Px() + ptkick_*charge;
	double py = vec.Py();
	double pz = sqrt(p*p - px*px - py*py);

	TLorentzVector kicked;
	kicked.SetXYZM(px,py,pz,vec.M());

	return kicked;
}
