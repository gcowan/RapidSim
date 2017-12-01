#ifndef RAPIDPID_H
#define RAPIDPID_H

#include <map>
#include <set>

#include "TH3.h"
#include "TString.h"

class RapidPID {
	public:
		RapidPID(TString name)
			: name_(name) {}

		~RapidPID();

		double getPID(unsigned int id, double p, double eta);
		void addPID(unsigned int id, TH3D* hist);

	private:
		void setupRangePEta(unsigned int id);
		void limitRangePEta(unsigned int id, double& p, double& eta);

		TString name_;

		std::map<unsigned int, TH3D*> pidHists_;
		std::map<unsigned int, std::map<unsigned int, TH1D*>* > cachedPIDHists_;

		std::map<unsigned int, double> maxP_;
		std::map<unsigned int, double> minEta_;
		std::map<unsigned int, double> maxEta_;

		std::set<unsigned int> suppressWarning_;
};

#endif
