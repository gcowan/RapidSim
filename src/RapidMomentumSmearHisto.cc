#include "RapidMomentumSmearHisto.h"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

RapidMomentumSmearHisto::~RapidMomentumSmearHisto() {
	while(!histos_.empty()) {
		delete histos_[histos_.size()-1];
		histos_.pop_back();
	}
}

TLorentzVector RapidMomentumSmearHisto::smearMomentum(TLorentzVector p) {


	double kp, kpT, kptx, kpty, norm, smear;
	kp = p.P();
	kpT = p.Pt();
	kptx = p.Px()/p.Pz();
	kpty = p.Py()/p.Pz();
	unsigned int iHist=0;
	while( true ) {
		if( iHist==thresholds_.size()-1 ) break;
		// if( kp < thresholds_[iHist+1] ) break;
		if( kpT < thresholds_[iHist+1] ) break;
		++iHist;
	}

	// std::cout << "SMEARING" << " " << iHist << " " << kp<< std::endl;

	smear = histos_[iHist]->GetRandom()*kp;
	//smear = 1.0*ran.Gaus(0,1)*dGraph->Eval(1000*kp)*kp;
	kp += smear;

	// smear the slopes
	double slope_smear = TMath::Sqrt(TMath::Power(6.2e-5,2) + TMath::Power(2.1e-3/kp,2)); //TODO
	kptx += slope_smear*gRandom->Gaus(1,0);
	kpty += slope_smear*gRandom->Gaus(1,0);
	norm = sqrt(1 + kptx*kptx + kpty*kpty);
	if(p.Pz()<0) norm = -norm;

	p.SetXYZM( kptx*kp/norm, kpty*kp/norm, kp/norm, p.M() );
	return p;

}

void RapidMomentumSmearHisto::init(std::vector<double> thresholds, std::vector<TH1*> histos) {
	if(thresholds.size() < histos.size()) {
		std::cout << "WARNING in RapidMomentumSmearHisto::init : too many histograms provided. Number of histograms should match number of thresholds." << std::endl;
		std::cout << "                                      excess histograms ignored." << std::endl;

		while(thresholds.size() < histos.size()) {
			histos.pop_back();
		}
	} else if(thresholds.size() > histos.size()) {
		std::cout << "WARNING in RapidMomentumSmearHisto::init : too few histograms provided. Number of histograms should match number of thresholds." << std::endl;
		std::cout << "                                      excess thresholds ignored." << std::endl;

		while(thresholds.size() > histos.size()) {
			thresholds.pop_back();
		}
	}

	thresholds_ = thresholds;
	histos_ = histos;
}
