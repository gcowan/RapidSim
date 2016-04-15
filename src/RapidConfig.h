#ifndef RAPIDCONFIGPARSER_H
#define RAPIDCONFIGPARSER_H

#include <map>
#include <vector>

#include "TH1F.h"
#include "TString.h"

class RapidMomentumSmear;
class RapidParam;
class RapidParticle;

class RapidConfig {
	public:
		RapidConfig()
			: fileName_(""), accRejHisto_(0), accRejParameter_(0)
		{}

		~RapidConfig();

		void load(TString fileName);

		const std::vector<RapidParticle*>& particles() { return parts_; }
		const std::vector<RapidParam*>& parameters() { return params_; }
		
		TH1F* acceptRejectHistogram() { return accRejHisto_; }
		RapidParam* acceptRejectParameter() { return accRejParameter_; }

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
};

#endif
