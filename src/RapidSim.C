#include <iostream>
#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "functions.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooRelBreitWigner.h"

#include "RapidDecay.h"
#include "RapidParticleData.h"

void rapidSim(const std::string mode, const int nEvtToGen, const std::string path, bool saveTree=false) {

	// load the fonll stuff
	TFile* fonll = new TFile((path + "/fonll/fonll.root").c_str());
	TH1F* ptHisto = (TH1F*) fonll->Get("pthisto");
	TH1F* etaHisto = (TH1F*) fonll->Get("etahisto");

	TRandom3 ran;

	RapidParticleData* rpd = RapidParticleData::getInstance();
	rpd->loadData("../config/particles.dat");

	RapidDecay myDecayObject(mode);
	myDecayObject.setRandomGenerator(ran);
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

	const std::string mode = argv[1];
	const int number = atoi(argv[2]);
	const std::string path = argv[3];
	bool saveTree = false;

	if(argc>4) {
		saveTree = atoi(argv[4]);
	}

	rapidSim(mode, number, path, saveTree);
	return 0;
}
