#ifndef RAPIDDECAY_H
#define RAPIDDECAY_H

#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "RapidConfig.h"
#include "RapidHistWriter.h"
#include "RapidParticleData.h"

class RapidParticle;
class RapidParam;

class RapidDecay {
	public:
		RapidDecay(TString filename, bool saveTree)
			: maxgen_(1000),
			  ptHisto_(0), etaHisto_(0), accRejHisto_(0),
			  particleData_(RapidParticleData::getInstance())
			{loadDecay(filename, saveTree);}

		~RapidDecay() {}
		
		void setMaxGen(int mg) { maxgen_ = mg; }
		void loadParentKinematics(TH1F* ptHisto, TH1F* etaHisto);
		
		bool generate();

		void saveHistos();

	private:
		void loadDecay(TString filename, bool saveTree);

		void setAcceptRejectHist(TH1F* histo, RapidParam* param);
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

		//particle data lookup
		RapidParticleData* particleData_;

		//config loader
		RapidConfig config_;

		//histogram writer
		RapidHistWriter writer_;
};
#endif
