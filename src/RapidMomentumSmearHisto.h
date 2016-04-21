#ifndef RAPIDMOMENTUMSMEARHISTO_H
#define RAPIDMOMENTUMSMEARHISTO_H

#include <vector>

#include "TH1.h"

#include "RapidMomentumSmear.h"

class RapidMomentumSmearHisto : public RapidMomentumSmear {
	public:
		RapidMomentumSmearHisto(std::vector<double> thresholds, std::vector<TH1*> histos) { init(thresholds, histos); }

		~RapidMomentumSmearHisto();

		TLorentzVector smearMomentum(TLorentzVector p);

	private:
		void init(std::vector<double> thresholds, std::vector<TH1*> histos);

		std::vector<double> thresholds_;
		std::vector<TH1*> histos_;
};

#endif
