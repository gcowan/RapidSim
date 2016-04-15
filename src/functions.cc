#include "functions.h"

TLorentzVector toFourVector(const TVector3& vec, double m) {
    return TLorentzVector(vec, TMath::Sqrt(m*m + vec.Mag2()));
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

    return (TMath::Abs(xTracker) < xSizeTracker && TMath::Abs(xTracker) > xMinTracker 
            &&  TMath::Abs(yTracker) < ySizeTracker) && TMath::Abs(yTracker) > yMinTracker ;  
}

