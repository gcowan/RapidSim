#ifndef RAPIDHISTWRITER_H
#define RAPIDHISTWRITER_H

#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"

class RapidParam;
class RapidParticle;

class RapidHistWriter {
	public:
		RapidHistWriter(const std::vector<RapidParticle*>& parts, const std::vector<RapidParam*>& params, const std::vector<RapidParam*>& paramsStable, const std::vector<RapidParam*>& paramsDecaying, const std::vector<RapidParam*>& paramsTwoBody, const std::vector<RapidParam*>& paramsThreeBody, TString name, bool saveTree)
			: name_(name), parts_(parts), params_(params), paramsStable_(paramsStable), paramsDecaying_(paramsDecaying), paramsTwoBody_(paramsTwoBody), paramsThreeBody_(paramsThreeBody), treeFile_(0), tree_(0)
		{setup(saveTree);}

		~RapidHistWriter();

		void fill();
		void save();

		void setNEvent(int nevent) { nevent_ = nevent; }

	private:
		void setup(bool saveTree);

		void setupHistos();
		void setupSingleHypothesis(TString suffix="");
		void setupTree();

		unsigned int fillSingleHypothesis(unsigned int offset=0);

		TString name_;

		std::vector<RapidParticle*> parts_;

		//particles with alternative mass hypotheses
		std::vector<RapidParticle*> altHypothesisParts_;

		//particle specific parameters
		std::vector<RapidParam*> params_;

		//default parameter sets
		std::vector<RapidParam*> paramsStable_;
		std::vector<RapidParam*> paramsDecaying_;
		std::vector<RapidParam*> paramsTwoBody_;
		std::vector<RapidParam*> paramsThreeBody_;

		//histograms to store parameters in
		std::vector<TH1F*> histos_;

		//tree to store parameters in
		TFile* treeFile_;
		TTree* tree_;
		int nevent_;
		std::vector<double> vars_;
};

#endif
