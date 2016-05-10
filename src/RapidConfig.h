#ifndef RAPIDCONFIGPARSER_H
#define RAPIDCONFIGPARSER_H

#include <map>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"

#include "RapidAcceptance.h"

class RapidCut;
class RapidDecay;
class RapidHistWriter;
class RapidMomentumSmear;
class RapidParam;
class RapidParticle;

class RapidConfig {
	public:
		RapidConfig()
			: fileName_(""), accRejHisto_(0),
			  accRejParameterX_(0), accRejParameterY_(0),
			  acceptanceType_(RapidAcceptance::ANY),
			  ppEnergy_(8.), motherFlavour_("b"),
			  ptHisto_(0), etaHisto_(0), maxgen_(1000),
			  decay_(0), acceptance_(0), writer_(0)
		{}

		~RapidConfig();

		bool load(TString fileName);

		RapidDecay* getDecay();
		RapidAcceptance* getAcceptance();
		RapidHistWriter* getWriter(bool saveTree=false);

	private:
		bool loadDecay();
		bool loadConfig();
		void writeConfig();

		bool configParticle(unsigned int part, TString command, TString value);
		bool configGlobal(TString command, TString value);
		RapidParam* loadParam(TString paramStr);
		RapidCut* loadCut(TString cutStr);

		RapidParam* findParam(TString name);

		void setSmearing(unsigned int particle, TString category);
		bool loadSmearing(TString category);

		bool loadAcceptRejectHist(TString histFile, TString histName, RapidParam* paramX, RapidParam* paramY);
		bool loadParentKinematics();

		bool check1D(TH1* hist) { return (dynamic_cast<TH1F*>(hist) || dynamic_cast<TH1D*>(hist)); }
		bool check2D(TH1* hist) { return (dynamic_cast<TH2F*>(hist) || dynamic_cast<TH2D*>(hist)); }

		TString fileName_;

		std::vector<RapidParticle*> parts_;
		std::vector<RapidParam*> params_;

		std::vector<RapidCut*> cuts_;

		//mometum smearing lookup for each smearing category
		std::map<TString, RapidMomentumSmear*> momSmearCategories_;

		//accept reject hist to sculpt kinematics
		TH1* accRejHisto_;
		RapidParam* accRejParameterX_;
		RapidParam* accRejParameterY_;

		//type of geometric acceptance to apply
		RapidAcceptance::AcceptanceType acceptanceType_;

		//parameters to determine parent kinematics
		double ppEnergy_;
		TString motherFlavour_;

		//parent kinematic distributions
		TH1* ptHisto_;
		TH1* etaHisto_;

		//max attempts to generate
		double maxgen_;

		RapidDecay* decay_;
		RapidAcceptance* acceptance_;
		RapidHistWriter* writer_;
};

#endif
