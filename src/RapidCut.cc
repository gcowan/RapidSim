#include "RapidCut.h"

#include "RapidParam.h"

const double RapidCut::NOLIMIT = -99999.;

bool RapidCut::passCut() {
	double val = param_->eval();

	bool inRange = (val>min_ || min_==NOLIMIT) && (val<max_ || max_==NOLIMIT);

	return (inRange!=veto_);

}

TString RapidCut::name() {
	TString name("");
	if(veto_) {
		TString vetoStr("");
		vetoStr.Form(" veto(%.3f:%.3f)", min_, max_);
		name += param_->name();
		name += vetoStr;
		//name += " veto(";
		//name += min_;
		//name += ":";
		//name += max_;
		//name += ")";
	} else {
		TString minStr(""), maxStr("");
		if(min_!=NOLIMIT) {
			minStr.Form("%.3f < ",min_);
			//name += min_;
			//name += " < ";
		}
		name += minStr;
		name += param_->name();
		if(max_!=NOLIMIT) {
			maxStr.Form(" < %.3f",max_);
			//name += " < ";
			//name += max_;
		}
		name += maxStr;
	}

	return name;
}
