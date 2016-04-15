#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

#include "RapidDecay.h"
#include "RapidParticleData.h"

void rapidSim(const TString mode, const int nEvtToGen, const TString path, bool saveTree=false) {

	// load the fonll stuff
	TFile* fonll = new TFile(path + "/fonll/fonll.root");
	TH1F* ptHisto = (TH1F*) fonll->Get("pthisto");
	TH1F* etaHisto = (TH1F*) fonll->Get("etahisto");

	RapidParticleData* rpd = RapidParticleData::getInstance();
	rpd->loadData("../config/particles.dat");

	RapidDecay myDecayObject(mode);
	myDecayObject.loadParentKinematics(ptHisto,etaHisto);

	if(saveTree) myDecayObject.saveTree("myTree.root");


	int ngenerated = 0; int nselected = 0;
	for (Int_t n=0; n<nEvtToGen; ++n) {
		if (!myDecayObject.generate()) continue;
		++ngenerated;
		++nselected;
	} //event loop

	myDecayObject.saveHistos("myHistos.root");

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
