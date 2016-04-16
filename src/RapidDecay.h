#ifndef RAPIDDECAY_H
#define RAPIDDECAY_H

#include <vector>

#include "TH1F.h"
#include "TString.h"

class RapidParticle;
class RapidParam;

class RapidDecay {
	public:
		RapidDecay(const std::vector<RapidParticle*>& parts)
			: parts_(parts), maxgen_(1000),
			  ptHisto_(0), etaHisto_(0),
			  accRejHisto_(0), accRejParameter_(0)
			{setup();}

		~RapidDecay() {}
		
		void setMaxGen(int mg) { maxgen_ = mg; }
		void setParentKinematics(TH1F* ptHisto, TH1F* etaHisto);
		void setAcceptRejectHist(TH1F* histo, RapidParam* param);
		
		bool generate();

	private:
		void setup();

		bool loadSmearing(TString category);
		void setupMasses();

		bool runAcceptReject();
		TH1F* generateAccRejDenominator();

		void floatMasses();
		void genParent();
		bool genDecay(bool acceptAny=false);
		bool genDecayAccRej();
		void smearMomenta();

		//the particles
		std::vector<RapidParticle*> parts_;
		
		//max number of attempts to generate an event
		int maxgen_;

		//parent kinematics
		TH1F* ptHisto_;
		TH1F* etaHisto_;

		//accept reject hist to sculpt kinematics
		TH1F* accRejHisto_;
		RapidParam* accRejParameter_;

};
#endif
