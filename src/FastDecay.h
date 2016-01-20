#ifndef DECAY_H
#define DECAY_H

#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "RooDataSet.h"

#include <iostream>
#include <vector>
#include <map>

#include "functions.h"

class FastDecay {
	public:
		enum ParamType {
			M,        //Mass
			M2,       //Mass squared
			MT,       //Transverse mass
			E,        //Energy
			ET,       //Transverse energy
			P,        //3-momentum
			PX,       //X-momentum
			PY,       //Y-momentum
			PZ,       //Z-momentum
			PT,       //Transverse momentum
			ETA,      //Pseudorapidity
			PHI,      //Azimuthal angle
			RAPIDITY, //Rapidity
			GAMMA,    //Relitivistic gamma
			BETA      //velocity

		};

		FastDecay(TString filename)
			: tree(0), varsPerPart(0),
			  rand(0), maxgen(1000),
			  ptHisto(0), etaHisto(0), smearGraph(0), accRejHisto(0)
			{loadDecay(filename);}

		void setRandomGenerator(TRandom& r) { rand = r; }
		void setMaxGen(int mg) { maxgen = mg; }
		void loadParentKinematics(TH1F* pt, TH1F* eta);
		void loadSmearGraph(TGraphErrors* sg) { smearGraph = sg;}
		void setAcceptRejectHist(TH1F* hist, ParamType type, std::vector<int> particles);
		void addCustomParameter(TString name, ParamType type, std::vector<int> particles, bool truth=false, double min=0., double max=10.);
		
		bool generate();
		TLorentzVector getP(unsigned int i);
		TLorentzVector getPSmeared(unsigned int i);
		void saveHistos(TString fname);
		void saveTree(TString fname);

	protected:
		struct CustomParameter {
			TString name;
			ParamType type;
			std::vector<int> particles;
			bool truth;
			double minVal;
			double maxVal;
		};
	private:
		void loadDecay(TString filename);
		void setupMasses();
		TString getUniqName(TString base);
		void setupHistos();
		void setupTree();

		double evalCustomParam(int i);
		double evalCustomParam(CustomParameter param);

		bool runAcceptReject();
		TH1F* generateAccRejDenominator();

		void floatMasses();
		void genParent();
		bool genDecay();
		void smearMomenta();
		void fillHistos();
		void fillTree();

		void setupRhoMass();
		void setupKstMass();
		void setupPhiMass();
		void setupChic0Mass();
		void setupChic1Mass();
		void setupChic2Mass();

		std::vector<int> parts;
		std::vector<TString> names;
		std::vector<int> mother;
		std::vector<int> nDaug;
		std::vector<int> firstDaug;
		std::vector<double> m;
		std::vector<TLorentzVector> p;
		std::vector<TLorentzVector> pSmeared;
		std::vector<TH1F*> histos;

		//register of particle names that have been used
		std::set<TString> usedNames;

		//tree to store parameters in
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

		//momentum smearing
    		TGraphErrors* smearGraph;

		//accept reject hist to sculpt kinematics
		TH1F* accRejHisto;
		CustomParameter accRejParameter;

		//custom parameters added to the histograms and the tree
		std::vector<CustomParameter> customParams;

		//datasets to sample resonance masses from
		//loaded on-demand
		std::map<int, RooDataSet*> massdata;
		std::map<int, double> minmass;
		std::map<int, double> maxmass;
};
#endif
