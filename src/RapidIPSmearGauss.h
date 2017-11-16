#ifndef RAPIDIPSMEARGAUSS_H
#define RAPIDIPSMEARGAUSS_H

#include "TGraphErrors.h"

#include "RapidIPSmear.h"

class RapidIPSmearGauss : public RapidIPSmear {
	public:
		RapidIPSmearGauss(double intercept, double slope) : intercept_(intercept),slope_(slope) {}

		~RapidIPSmearGauss() {}

		std::pair<double,double> smearIP(double ip, double pt);

	private:
        double intercept_,slope_;

};

#endif
