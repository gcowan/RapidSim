#include "RapidMomentumSmearBrems.h"
#include "TMath.h"
#include "TRandom.h"

RapidMomentumSmearBrems::~RapidMomentumSmearBrems() {
  if (lhcb_) delete lhcb_;
  if(graph_) delete graph_;
}

TLorentzVector RapidMomentumSmearBrems::smearMomentum(TLorentzVector p) {

    Event event;
    event.push_back(Particle(p.Px(), p.Py(), p.Pz(), sqrt(p.P()*p.P() + 5.110e-04*5.110e-04), 11));
    event[0].x[0] = 0.46;
    event[0].x[1] = -0.013;
    event[0].x[2] = 0.94; 
    lhcb_->transport(event);

    double kp, kptx, kpty, norm, smear;
    kp = sqrt(event[0][0]*event[0][0]+event[0][1]*event[0][1]+event[0][2]*event[0][2]);
    kptx = event[0][0]/event[0][2];
    kpty = event[0][1]/event[0][2];
    smear = 1.0*gRandom->Gaus(0,1)*graph_->Eval(1000*kp)*kp;
    kp += smear;

    double slope_smear = TMath::Sqrt(TMath::Power(6.2e-5,2) + TMath::Power(2.1e-3/kp,2)); //TODO
    kptx += slope_smear*gRandom->Gaus(1,0);
    kpty += slope_smear*gRandom->Gaus(1,0);
    norm = sqrt(1 + kptx*kptx + kpty*kpty);
    if(event[0][2]<0) norm = -norm;

    p.SetXYZM( kptx*kp/norm, kpty*kp/norm, kp/norm, p.M() );
    return p;

}

