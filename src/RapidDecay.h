#ifndef RAPIDDECAY_H
#define RAPIDDECAY_H

#include <vector>
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"

#include "RapidVertex.h"

class RapidParticle;
class RapidParam;
class RapidExternalGenerator;

class RapidDecay {
	public:
		RapidDecay(const std::vector<RapidParticle*>& parts)
			: parts_(parts), maxgen_(1000),
			  ptHisto_(0), etaHisto_(0),
			  pvHisto_(0),
			  accRejHisto_(0), accRejParameterX_(0), accRejParameterY_(0),
			  suppressKinematicWarning_(false), suppressAttemptsWarning_(false),
			  external_(0)
			{setup();}

		~RapidDecay() {}

		void setMaxGen(int mg) { maxgen_ = mg; }
		void setParentKinematics(TH1* ptHisto, TH1* etaHisto);
		void setParentKinematics2d(TH2* momentaHisto2D);
		void setPVntracks(TH1* pvHisto);
		void setAcceptRejectHist(TH1* histo, RapidParam* param);
		void setAcceptRejectHist(TH1* histo, RapidParam* paramX, RapidParam* paramY);
		void setExternal(RapidExternalGenerator* external);

		bool checkDecay();
		bool generate(bool genpar=true);

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
		void calcIPs();
		double getParticleIP(ROOT::Math::XYZPoint, ROOT::Math::XYZPoint, TLorentzVector);

		//the particles
		std::vector<RapidParticle*> parts_;

		//pileup vertices
		std::vector<RapidVertex> pileuppvs_;

		//max number of attempts to generate an event
		int maxgen_;

		//parent kinematics
		TH1* ptHisto_;
		TH1* etaHisto_;

		TH2* momentaHisto2D_;

		//PVNTRACKS
		TH1* pvHisto_;

		//accept reject hist to sculpt kinematics
		TH1* accRejHisto_;
		RapidParam* accRejParameterX_;
		RapidParam* accRejParameterY_;

		//TGenPhaseSpace object to perform the decays
		TGenPhaseSpace decay_;

		//flags to suppress generation warnings
		bool suppressKinematicWarning_;
		bool suppressAttemptsWarning_;

		//external decay generator wrapper
		RapidExternalGenerator* external_;

};
#endif
