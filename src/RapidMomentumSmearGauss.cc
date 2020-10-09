#include "RapidMomentumSmearGauss.h"

#include "TMath.h"
#include "TRandom.h"

RapidMomentumSmearGauss::~RapidMomentumSmearGauss() {
	if(graph_) delete graph_;
}

TLorentzVector RapidMomentumSmearGauss::smearMomentum(TLorentzVector p) {

	double kp, kptx, kpty, norm, smear;
	kp = p.P();
	kptx = p.Px()/p.Pz();
	kpty = p.Py()/p.Pz();
	smear = 1.0*gRandom->Gaus(0,1)*graph_->Eval(1000*kp)*kp;
	kp += smear;

	// smear the slopes
	double slope_smear = TMath::Sqrt(TMath::Power(6.2e-5,2) + TMath::Power(2.1e-3/kp,2)); //TODO
	kptx += slope_smear*gRandom->Gaus(0,1);
	kpty += slope_smear*gRandom->Gaus(0,1);
	norm = sqrt(1 + kptx*kptx + kpty*kpty);
	if(p.Pz()<0) norm = -norm;

	p.SetXYZM( kptx*kp/norm, kpty*kp/norm, kp/norm, p.M() );
	return p;
}

