#ifndef RAPIDCUT_H
#define RAPIDCUT_H

#include "TString.h"

class RapidParam;

class RapidCut {
	public:
		RapidCut(RapidParam* param, double min=NOLIMIT, double max=NOLIMIT, bool veto=false)
			: param_(param), min_(min), max_(max), veto_(veto) {}
		
		TString name();

		bool passCut();

		static const double NOLIMIT;

	private:
		RapidParam* param_;
		double min_;
		double max_;
		bool veto_;
};

#endif
