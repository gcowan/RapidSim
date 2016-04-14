#include "functions.h"

double resSlope(double p) {
    // std::cout << p << std::endl;
    return(sqrt(pow(6.2e-5,2) + pow(2.1e-3/p,2)));  
}

TLorentzVector toFourVector(const TVector3& vec, double m) {
    return TLorentzVector(vec, TMath::Sqrt(m*m + vec.Mag2()));
}

TLorentzVector reassignMass(const TLorentzVector& vec, double mass){
    const TVector3 threevec = vec.Vect();
    return toFourVector(threevec, mass);
}

TVector3 smearedVec(double plus_x, double plus_y, double plus_z, TGraphErrors* dGraph, TRandom& ran) {
    double kp, kptx, kpty, norm, smear;
    kp = sqrt(plus_x*plus_x + plus_y*plus_y + plus_z*plus_z );
    kptx = plus_x/plus_z;
    kpty = plus_y/plus_z;
    norm = sqrt(1 + kptx*kptx + kpty*kpty);
    smear = 1.0*ran.Gaus(0,1)*dGraph->Eval(1000*kp)*kp;
    kp += smear;

    // smear the slopes
    double slope_smear = resSlope(kp); 
    kptx += slope_smear*ran.Gaus(1,0);
    kpty += slope_smear*ran.Gaus(1,0);
    norm = sqrt(1 + kptx*kptx + kpty*kpty);

    plus_x = kptx*kp/norm;
    plus_y = kpty*kp/norm;
    plus_z = kp/norm;
    return TVector3(plus_x,plus_y, plus_z);

}

TLorentzVector smearedVec(TLorentzVector& vec, TGraphErrors* dGraph, TRandom& ran){

    TVector3 threeVec = smearedVec(vec.Px(), vec.Py(), vec.Pz(),dGraph,ran);
    return toFourVector(threeVec, vec.M());
}

bool inAcceptance(TLorentzVector& vec){

    if (TMath::Abs(vec.Px()/vec.Pz()) > 0.3) return false;
    if (TMath::Abs(vec.Py()/vec.Pz()) > 0.25) return false;  
    if (sqrt(pow(vec.Px()/vec.Pz(),2) + pow(vec.Py()/vec.Pz(),2)) <0.01) return false;

    return true;
}

TLorentzVector magnetKick(TLorentzVector& vec, double charge){

    TVector3 threeVec  =  vec.Vect();

    // extrapolate to magnet centre  

    // kick 
    double p = vec.P();
    double px =  vec.Px() + ptkick*charge ;
    double py = vec.Py();
    double pz = sqrt(p*p - px*px - py*py);

    TVector3 threevec = TVector3(px,py,pz) ;

    return toFourVector(threevec, vec.M()); 

}

bool inDownstream(TLorentzVector& vec, int charge){

    TLorentzVector newvec = magnetKick(vec,charge);

    // position at magnet centre
    double xMag = zC*vec.Px()/vec.Pz(); 
    double yMag = zC*vec.Py()/vec.Pz();  

    double xTracker = xMag +  (newvec.Px()*(zTracker - zC)/newvec.Pz()); 
    double yTracker = yMag +  (newvec.Py()*(zTracker - zC)/newvec.Pz());   

    //if (TMath::Abs(xTracker) > xSizeTracker) std::cout << "out of tracker " << std::endl;

    return (TMath::Abs(xTracker) < xSizeTracker && TMath::Abs(xTracker) > xMinTracker 
            &&  TMath::Abs(yTracker) < ySizeTracker) && TMath::Abs(yTracker) > yMinTracker ;  
}

