#include "RapidMomentumSmearEnergyGauss.h"

#include "TMath.h"
#include "TRandom.h"
#include <iostream>

TLorentzVector RapidMomentumSmearEnergyGauss::smearMomentum(TLorentzVector p) {
    double energy = p.E();
    double first = stochastic_/TMath::Sqrt(energy);
    first *= first;
    double second = constant_*constant_;
    double res = TMath::Sqrt(first + second)*energy;
    double smearedEnergy = gRandom->Gaus(0,res) + energy;
    double norm = smearedEnergy/energy;
	p.SetPxPyPzE( p.Px()*norm, p.Py()*norm, p.Pz()*norm, smearedEnergy );
	return p;
}

