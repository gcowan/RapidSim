#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

#include "RapidAcceptance.h"
#include "RapidConfig.h"
#include "RapidDecay.h"
#include "RapidHistWriter.h"

void rapidSim(const TString mode, const int nEvtToGen, const TString path, bool saveTree=false) {

	// load the fonll stuff
	TFile* fonll = new TFile(path + "/fonll/fonll.root");
	TH1F* ptHisto = (TH1F*) fonll->Get("pthisto");
	TH1F* etaHisto = (TH1F*) fonll->Get("etahisto");

	RapidConfig config;
	config.load(mode);

	RapidDecay* decay = config.getDecay();
	decay->setParentKinematics(ptHisto,etaHisto);

	RapidAcceptance* acceptance = config.getAcceptance();

	RapidHistWriter* writer = config.getWriter(saveTree);

	int ngenerated = 0; int nselected = 0;
	for (Int_t n=0; n<nEvtToGen; ++n) {
		if (!decay->generate()) continue;
		++ngenerated;

		if(!acceptance->inAcceptance()) continue;
		++nselected;

		writer->fill();
	}

	writer->save();

	std::cout << "Generated " << ngenerated << std::endl;
	std::cout << "Selected " << nselected << std::endl;

}

int main(int argc, char * argv[])
{
	if (argc < 4) {
		printf("Usage: %s mode numberToGenerate pathToFiles [saveTree=false]\n", argv[0]);
		return 1;
	}

	const TString mode = argv[1];
	const int number = atoi(argv[2]);
	const TString path = argv[3];
	bool saveTree = false;

	if(argc>4) {
		saveTree = atoi(argv[4]);
	}

	rapidSim(mode, number, path, saveTree);
	return 0;
}
