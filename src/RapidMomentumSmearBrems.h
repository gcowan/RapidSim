#ifndef RAPIDMOMENTUMSMEARBREMS_H
#define RAPIDMOMENTUMSMEARBREMS_H

#include "TGraphErrors.h"
#include "fastsim.h"
#include "RapidMomentumSmear.h"

class RapidMomentumSmearBrems : public RapidMomentumSmear {
	public:
		RapidMomentumSmearBrems(Detector* lhcb,TGraphErrors* graph) : lhcb_(lhcb), graph_(graph) {}

		~RapidMomentumSmearBrems();

		TLorentzVector smearMomentum(TLorentzVector p);
    private:
      Detector* lhcb_;
      TGraphErrors* graph_;

};

#endif
