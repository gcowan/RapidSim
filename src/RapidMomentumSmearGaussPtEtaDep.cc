#include "RapidMomentumSmearGaussPtEtaDep.h"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

RapidMomentumSmearGaussPtEtaDep::~RapidMomentumSmearGaussPtEtaDep() {
	if(hist_) delete hist_;
}

TLorentzVector RapidMomentumSmearGaussPtEtaDep::smearMomentum(TLorentzVector p) {

	double smear, pt, eta;
	pt = p.Pt();
	eta = p.Eta();

	int bin = hist_->FindBin(pt,eta);

	smear = 1.0*gRandom->Gaus()*hist_->GetBinContent(bin);

	p.SetXYZM( p.Px()*(1+smear), p.Py()*(1+smear), p.Pz(), p.M() );
	return p;

}
