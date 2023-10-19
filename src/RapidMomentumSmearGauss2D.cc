#include "RapidMomentumSmearGauss2D.h"

#include "TMath.h"
#include "TRandom.h"

RapidMomentumSmearGauss2D::~RapidMomentumSmearGauss2D() {
	if(pRes_)  delete pRes_;
	if(txRes_) delete txRes_;
	if(tyRes_) delete tyRes_;
}
TLorentzVector RapidMomentumSmearGauss2D::smearMomentum(TLorentzVector p) {
	double kp, kptx, kpty, norm, smear , kpt, eta;
	kp   = p.P();
	kpt  = TMath::Sqrt( p.Px()*p.Px() + p.Py() * p.Py());
	kptx = p.Px()/p.Pz();
	kpty = p.Py()/p.Pz();
	// p-resolution parameterisation is in x = log10(pt) , y = eta 
	eta = p.PseudoRapidity();
	auto getZFromXY = [&]( TProfile2D * profile, double inputX, double inputY){
		// std::cout<< "X, Y input = "<< inputX << ","<< inputY << std::endl;
		double use_inputX = inputX;
		double use_inputY = inputY;
		auto binIDX = profile->FindFixBin( inputX, inputY);	
		if(profile->IsBinOverflow(binIDX) || profile->IsBinUnderflow(binIDX)){
			double Res_xLow = profile->GetXaxis()->GetXmin();
			double Res_xHigh= profile->GetXaxis()->GetXmax();
			double Res_yLow = profile->GetYaxis()->GetXmin();
			double Res_yHigh= profile->GetYaxis()->GetXmax();
			// std::cout<<"x Limits = ["<< Res_xLow <<","<< Res_xHigh << "]" << std::endl;
			// std::cout<<"y Limits = ["<< Res_yLow <<","<< Res_yHigh << "]" << std::endl;
			if( inputX < Res_xLow)  use_inputX = Res_xLow+0.00001;
			if( inputX > Res_xHigh) use_inputX = Res_xHigh-0.00001;
			if( inputY < Res_yLow)  use_inputY = Res_yLow+0.00001;
			if( inputY > Res_yHigh) use_inputY = Res_yHigh-0.00001;
			// std::cout<< "X, Y input (new) = "<< use_inputX << ","<< use_inputY << std::endl;			
			return profile->Interpolate( use_inputX , use_inputY) ; 
		}
		return profile->Interpolate( use_inputX , use_inputY) ; 
	};
	//Momentum smearing 
	double resolution =  (getZFromXY( pRes_, TMath::Log10( kpt*1000.), eta)/100. )* kp ;
	smear = 1.0*gRandom->Gaus(0,1) * resolution;
	kp+= smear;
	

	// Tx Slope smearing /Momentum smearing 
	// The histogram is done in bins of 1./log10(pt MeV) and eta
	double resolutionTX = 0.001*getZFromXY( txRes_, 1./TMath::Log10(kpt*1000.), eta); 
	double resolutionTY = 0.001*getZFromXY( tyRes_, 1./TMath::Log10(kpt*1000.), eta); 

	//maps are in mrad uncertainty, we need radians
	double slope_smearTX =  1.0*gRandom->Gaus(0,1) * resolutionTX;
	double slope_smearTY =  1.0*gRandom->Gaus(0,1) * resolutionTY;
	kptx += slope_smearTX;
	kpty += slope_smearTY;
	norm = sqrt(1 + kptx*kptx + kpty*kpty);
	if(p.Pz()<0) norm = -norm;
	p.SetXYZM( kptx*kp/norm, kpty*kp/norm, kp/norm, p.M() );
	return p;
}