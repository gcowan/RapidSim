#include "RapidVertex.h"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

ROOT::Math::XYZPoint RapidVertex::getVertex(bool truth) {
	if(truth) return vertexTrue_;
	else return vertexSmeared_;
}

void RapidVertex::setXYZ(double x, double y, double z) {
	vertexTrue_ = ROOT::Math::XYZPoint(x,y,z);
	smearVertex();
}

void RapidVertex::smearVertex() {
	// Obviously at the moment we are just using the same smearing for PV and SV.
	// Need to sample from some nPVTracks distribution
	// units are in mm
	double xS = 0.010817 + 0.03784*TMath::Exp(-0.0815*nPVTracks_);
	double yS = 0.010817 + 0.03784*TMath::Exp(-0.0815*nPVTracks_);
	double zS = 0.04252  + 0.2235 *TMath::Exp(-0.0814*nPVTracks_);
	vertexSmeared_ = ROOT::Math::XYZPoint( \
			vertexTrue_.X() + gRandom->Gaus(0,xS)*1000., \
			vertexTrue_.Y() + gRandom->Gaus(0,yS)*1000., \
			vertexTrue_.Z() + gRandom->Gaus(0,zS)*1000.);
}

