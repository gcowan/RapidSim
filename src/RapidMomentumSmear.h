#ifndef RAPIDMOMENTUMSMEAR_H
#define RAPIDMOMENTUMSMEAR_H

#include "TLorentzVector.h"

class RapidMomentumSmear {
	public:
		virtual TLorentzVector smearMomentum(TLorentzVector p)=0;
};

#endif
