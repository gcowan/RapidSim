#ifndef RAPIDBEAMDATA_H
#define RAPIDBEAMDATA_H

#include <map>
#include <set>

#include "TString.h"

#include "RooRealVar.h"

class RapidBeamData {
	public:

		static RapidBeamData* getInstance();

		void loadData(TString file);

		int getPileup() { return pileup_; }
		double getSigmaXY() { return sigmaxy_; }
		double getSigmaZ() { return sigmaz_; }

		
	private:

		static RapidBeamData* instance_;

        RapidBeamData()
		: pileup_(0), sigmaxy_(0.), sigmaz_(0.) {}

        ~RapidBeamData() {}

		//copy constructor and copy assignment operator not implemented
		RapidBeamData( const RapidBeamData& other );
		RapidBeamData& operator=( const RapidBeamData& other );

		int pileup_;
        double sigmaxy_;
        double sigmaz_;


};

#endif
