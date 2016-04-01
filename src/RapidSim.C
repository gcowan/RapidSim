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

void rapidSim(const std::string mode, const int nEvtToGen, const std::string path) {

    // load the fonll stuff
    TFile* fonll = new TFile((path + "/fonll/fonll.root").c_str());
    TH1F* ptHisto = (TH1F*) fonll->Get("pthisto"); 
    TH1F* etaHisto = (TH1F*) fonll->Get("etahisto"); 

    TRandom3 ran;

    RapidDecay myDecayObject(mode);
    myDecayObject.setRandomGenerator(ran);
    myDecayObject.loadParentKinematics(ptHisto,etaHisto);

//    std::vector<int> pars;
//    pars.push_back(2);
//    pars.push_back(3);
//    myDecayObject.addCustomParameter("qSqTRUE", RapidDecay::M2, pars, true, 0., 25.);
//    myDecayObject.addCustomParameter("qSq", RapidDecay::M2, pars, false, 0., 25.);
//    TH1F* arHist = new TH1F("arHist","",250,0.,25.);
//    for(int i=0; i<150; ++i) {
//	    arHist->SetBinContent(i+1, 1.);
//    }
//    myDecayObject.setAcceptRejectHist(arHist, RapidDecay::M2, pars);
    myDecayObject.saveTree("myTree.root");


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
    if (argc != 4) {
        printf("Usage: %s mode numberToGenerate pathToFiles\n", argv[0]);
        return 1;
    }
    const std::string mode = argv[1];
    const int number = atoi(argv[2]);
    const std::string path = argv[3];
    rapidSim(mode, number, path);
    return 0;
}
