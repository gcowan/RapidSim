#ifndef DECAY_H
#define DECAY_H

#include <iostream>
#include <vector>
#include <map>

#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "RooDataSet.h"

#include "functions.h"
#include "RapidParam.h"
#include "RapidParticleData.h"

class RapidParticle;
class RapidMomentumSmear;

class RapidDecay {
	public:
		RapidDecay(TString filename)
			: treeFile(0), tree(0), varsPerPart(0),
			  rand(0), maxgen(1000),
			  ptHisto(0), etaHisto(0), accRejHisto(0),
			  particleData(RapidParticleData::getInstance())
			{loadDecay(filename);}

		~RapidDecay();
		
		void setRandomGenerator(TRandom& r) { rand = r; }
		void setMaxGen(int mg) { maxgen = mg; }
		void loadParentKinematics(TH1F* pt, TH1F* eta);
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
		std::vector<RapidParticle*> parts;
		
		//custom parameters
		std::vector<TH1F*> histos;

		//tree to store parameters in
		TFile* treeFile;
		TTree* tree;
		std::vector<double> treeVars;
		int varsPerPart;

		//random number generator
		TRandom rand;

		//max number of attempts to generate an event
		int maxgen;

		//parent kinematics
    		TH1F* ptHisto; 
    		TH1F* etaHisto;

		//mometum smearing lookup for each smearing category
		std::map<TString, RapidMomentumSmear*> momSmearCategories;

		//accept reject hist to sculpt kinematics
		TH1F* accRejHisto;
		RapidParam* accRejParameter;

		//custom parameters added to the histograms and the tree
		std::vector<RapidParam*> customParams;

		//particle data lookup
		RapidParticleData* particleData;

};
#endif
