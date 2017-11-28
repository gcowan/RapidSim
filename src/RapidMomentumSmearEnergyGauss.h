#ifndef RAPIDMOMENTUMSMEARENERGYGAUSS_H
#define RAPIDMOMENTUMSMEARENERGYGAUSS_H

#include "RapidMomentumSmear.h"

class RapidMomentumSmearEnergyGauss : public RapidMomentumSmear {
	public:
		RapidMomentumSmearEnergyGauss(double stochastic, double constant) : stochastic_(stochastic), constant_(constant) {}

		TLorentzVector smearMomentum(TLorentzVector p);

	private:
		double stochastic_;
		double constant_;
};

#endif
