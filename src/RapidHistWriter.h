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
		RapidHistWriter()
			: name_(""), treeFile_(0), tree_(0), varsPerPart_(0) {}

		~RapidHistWriter();

		void setup(std::vector<RapidParticle*> parts, std::vector<RapidParam*> params, TString name, bool saveTree);
		void fill();
		void save();

	private:
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
