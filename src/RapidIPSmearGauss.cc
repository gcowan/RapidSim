#include "RapidIPSmearGauss.h"

#include "TMath.h"
#include "TRandom.h"

std::pair<double,double> RapidIPSmearGauss::smearIP(double ip, double pt) {

	const double sigma_ = intercept_ + slope_/pt;
	const double smear_ = ip+gRandom->Gaus(0.,sigma_);
	return std::pair<double,double>(smear_,sigma_);

}

