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
		RapidHistWriter(const std::vector<RapidParticle*>& parts, const std::vector<RapidParam*>& params, TString name, bool saveTree)
			: name_(name), parts_(parts), params_(params), treeFile_(0), tree_(0), varsPerPart_(0)
		{setup(saveTree);}

		~RapidHistWriter();

		void fill();
		void save();

	private:
		void setup(bool saveTree);

		void setupHistos();
		void setupTree();

		void fillHistos();
		void fillTree();

		TString name_;

		std::vector<RapidParticle*> parts_;
		std::vector<RapidParam*> params_;

		//histograms to store parameters in
		std::vector<TH1F*> histos_;

		//tree to store parameters in
		TFile* treeFile_;
		TTree* tree_;
		std::vector<double> treeVars_;
		int varsPerPart_;
};

#endif
