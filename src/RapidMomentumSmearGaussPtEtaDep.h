#ifndef RAPIDMOMENTUMSMEARGAUSSPTETADEP_H
#define RAPIDMOMENTUMSMEARGAUSSPTETADEP_H

#include <vector>

#include "TH2.h"

#include "RapidMomentumSmear.h"

class RapidMomentumSmearGaussPtEtaDep : public RapidMomentumSmear {
	public:
		RapidMomentumSmearGaussPtEtaDep(TH2* hist) : hist_(hist) {}

		~RapidMomentumSmearGaussPtEtaDep();

		TLorentzVector smearMomentum(TLorentzVector p);

	private:
		TH2* hist_;
};

#endif
