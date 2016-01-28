#ifndef RAPIDMOMENTUMSMEAR_H
#define RAPIDMOMENTUMSMEAR_H

#include <vector>

#include "TH1F.h"
#include "TLorentzVector.h"

class RapidMomentumSmear {
	public:
		RapidMomentumSmear(std::vector<double> thresholds, std::vector<TH1F*> histos) { init(thresholds, histos); }

		TLorentzVector smearMomentum(TLorentzVector p);

	private:
		void init(std::vector<double> thresholds, std::vector<TH1F*> histos);

		std::vector<double> thresholds_;
		std::vector<TH1F*> histos_;
};

#endif
