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

	if (this->partInAcceptance(part) == false) return false;

	TLorentzVector vec = part->getP();
	double charge = part->charge();
	TLorentzVector newvec = magnetKick(vec,charge);

	// position at origin
	ROOT::Math::XYZPoint orig = part->getOriginVertex()->getVertex(false);

	// position at magnet centre
	double xMag = 1e-3*orig.X() + (zC_ - 1e-3*orig.Z())*vec.Px()/vec.Pz();
	double yMag = 1e-3*orig.Y() + (zC_ - 1e-3*orig.Z())*vec.Py()/vec.Pz();

	double xTracker = xMag + (newvec.Px()*(zTracker_ - zC_)/newvec.Pz());
	double yTracker = yMag + (newvec.Py()*(zTracker_ - zC_)/newvec.Pz());

	bool inT = TMath::Abs(xTracker) < xSizeTracker_ && TMath::Abs(yTracker) < ySizeTracker_;
	bool inB = TMath::Abs(xTracker) < xMinTracker_ && TMath::Abs(yTracker) < yMinTracker_;

	return inT && !inB;
}

TLorentzVector RapidAcceptanceLHCb::magnetKick(TLorentzVector& vec, double charge){

	// kick
	double p = vec.P();
	double ty = vec.Px()/vec.Pz();
	double kick =  2.2*ty*ty - 0.005*ty + 1.26;
	if (kick > 1.4) kick = 1.4;
	double px = vec.Px() + kick*charge;
	double py = vec.Py();
	double pz = sqrt(p*p - px*px - py*py);

	TLorentzVector kicked;
	kicked.SetXYZM(px,py,pz,vec.M());

	return kicked;
}
