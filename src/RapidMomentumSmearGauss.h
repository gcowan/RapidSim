#ifndef RAPIDMOMENTUMSMEARGAUSS_H
#define RAPIDMOMENTUMSMEARGAUSS_H

#include <vector>

#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include "RapidMomentumSmear.h"

class RapidMomentumSmearGauss : public RapidMomentumSmear {
	public:
		RapidMomentumSmearGauss(TGraphErrors* graph) : graph_(graph) {}

		TLorentzVector smearMomentum(TLorentzVector p);

	private:
    		TGraphErrors* graph_;

};

#endif
