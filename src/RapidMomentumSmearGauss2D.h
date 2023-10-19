#ifndef RAPIDMOMENTUMSMEARGAUSS2D_H
#define RAPIDMOMENTUMSMEARGAUSS2D_H

#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "RapidMomentumSmear.h"

class RapidMomentumSmearGauss2D : public RapidMomentumSmear {
	public:
		RapidMomentumSmearGauss2D(TProfile2D* pRes, TProfile2D* txRes, TProfile2D* tyRes) : pRes_(pRes), txRes_(txRes), tyRes_(tyRes) {}
		~RapidMomentumSmearGauss2D();
		TLorentzVector smearMomentum(TLorentzVector p);
	private:
		TProfile2D* pRes_;
		TProfile2D* txRes_;
		TProfile2D* tyRes_;
};

#endif
