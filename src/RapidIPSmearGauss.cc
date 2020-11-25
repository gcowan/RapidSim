#include "RapidIPSmearGauss.h"

#include "TMath.h"
#include "TRandom.h"

std::pair<double,double> RapidIPSmearGauss::smearIP(double ip, double pt) {

	const double sigma_ = gRandom->Gaus(0.,intercept_ + slope_/pt);
	const double smear_ = ip+sigma_;
	return std::pair<double,double>(std::fabs(smear_),std::fabs(sigma_));

}

