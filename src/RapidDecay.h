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

#include "RapidParam.h"
#include "RapidParticleData.h"

class RapidParticle;
class RapidMomentumSmear;

class RapidDecay {
	public:
		RapidDecay(TString filename)
			: treeFile_(0), tree_(0), varsPerPart_(0),
			  maxgen_(1000),
			  ptHisto_(0), etaHisto_(0), accRejHisto_(0),
			  particleData_(RapidParticleData::getInstance())
			{loadDecay(filename);}

		~RapidDecay();
		
		void setMaxGen(int mg) { maxgen_ = mg; }
		void loadParentKinematics(TH1F* ptHisto, TH1F* etaHisto);
		void setSmearing(unsigned int particle, TString category);
		void setAcceptRejectHist(TString histFile, TString histName, RapidParam* param);
		
		bool generate();
		TLorentzVector getP(unsigned int i);
		TLorentzVector getPSmeared(unsigned int i);
		void saveHistos(TString fname);
		void saveTree(TString fname);

	private:
		void loadDecay(TString filename);
		void loadConfig(TString filename);
		void writeConfig(TString filename);
		void configParticle(unsigned int part, TString command, TString value);
		void configGlobal(TString command, TString value);
		RapidParam* loadParam(TString paramStr);

		bool loadSmearing(TString category);
		void setupMasses();
		TString getUniqName(TString base);
		void setupHistos();
		void setupTree();

		bool runAcceptReject();
		TH1F* generateAccRejDenominator();

		void floatMasses();
		void genParent();
		bool genDecay(bool acceptAny=false);
		bool genDecayAccRej();
		void smearMomenta();
		void fillHistos();
		void fillTree();

		//the particles
		std::vector<RapidParticle*> parts_;
		
		//custom parameters
		std::vector<TH1F*> histos_;

		//tree to store parameters in
		TFile* treeFile_;
		TTree* tree_;
		std::vector<double> treeVars_;
		int varsPerPart_;

		//max number of attempts to generate an event
		int maxgen_;

		//parent kinematics
		TH1F* ptHisto_;
		TH1F* etaHisto_;

		//mometum smearing lookup for each smearing category
		std::map<TString, RapidMomentumSmear*> momSmearCategories_;

		//accept reject hist to sculpt kinematics
		TH1F* accRejHisto_;
		RapidParam* accRejParameter_;

		//custom parameters added to the histograms and the tree
		std::vector<RapidParam*> customParams_;

		//particle data lookup
		RapidParticleData* particleData_;

};
#endif
