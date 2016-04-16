#ifndef RAPIDCONFIGPARSER_H
#define RAPIDCONFIGPARSER_H

#include <map>
#include <vector>

#include "TH1F.h"
#include "TString.h"

#include "RapidAcceptance.h"

class RapidDecay;
class RapidHistWriter;
class RapidMomentumSmear;
class RapidParam;
class RapidParticle;

class RapidConfig {
	public:
		RapidConfig()
			: fileName_(""), accRejHisto_(0), accRejParameter_(0),
			  acceptanceType_(RapidAcceptance::ANY),
			  decay_(0), acceptance_(0), writer_(0)
		{}

		~RapidConfig();

		void load(TString fileName);

		RapidDecay* getDecay();
		RapidAcceptance* getAcceptance();
		RapidHistWriter* getWriter(bool saveTree=false);

	private:
		void loadDecay();
		void loadConfig();
		void writeConfig();

		void configParticle(unsigned int part, TString command, TString value);
		void configGlobal(TString command, TString value);
		RapidParam* loadParam(TString paramStr);

		void setSmearing(unsigned int particle, TString category);
		bool loadSmearing(TString category);

		void loadAcceptRejectHist(TString histFile, TString histName, RapidParam* param);

		TString fileName_;

		std::vector<RapidParticle*> parts_;
		std::vector<RapidParam*> params_;

		//mometum smearing lookup for each smearing category
		std::map<TString, RapidMomentumSmear*> momSmearCategories_;

		//accept reject hist to sculpt kinematics
		TH1F* accRejHisto_;
		RapidParam* accRejParameter_;

		//type of geometric acceptance to apply
		RapidAcceptance::AcceptanceType acceptanceType_;

		RapidDecay* decay_;
		RapidAcceptance* acceptance_;
		RapidHistWriter* writer_;
};

#endif
