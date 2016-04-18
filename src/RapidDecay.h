#ifndef RAPIDDECAY_H
#define RAPIDDECAY_H

#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"

class RapidParticle;
class RapidParam;

class RapidDecay {
	public:
		RapidDecay(const std::vector<RapidParticle*>& parts)
			: parts_(parts), maxgen_(1000),
			  ptHisto_(0), etaHisto_(0),
			  accRejHisto_(0), accRejParameterX_(0), accRejParameterY_(0)
			{setup();}

		~RapidDecay() {}
		
		void setMaxGen(int mg) { maxgen_ = mg; }
		void setParentKinematics(TH1* ptHisto, TH1* etaHisto);
		void setAcceptRejectHist(TH1* histo, RapidParam* param);
		void setAcceptRejectHist(TH1* histo, RapidParam* paramX, RapidParam* paramY);
		
		bool generate();

	private:
		void setup();

		bool loadSmearing(TString category);
		void setupMasses();

		bool runAcceptReject();
		bool runAcceptReject1D();
		bool runAcceptReject2D();
		TH1* generateAccRejDenominator1D();
		TH2* generateAccRejDenominator2D();

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
		TH1* ptHisto_;
		TH1* etaHisto_;

		//accept reject hist to sculpt kinematics
		TH1* accRejHisto_;
		RapidParam* accRejParameterX_;
		RapidParam* accRejParameterY_;

};
#endif
