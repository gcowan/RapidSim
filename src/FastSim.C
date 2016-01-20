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
//#include "haofei.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooRelBreitWigner.h"

#include "FastDecay.h"

void PhaseHaofei(const int mode, const int nEvtToGen, const std::string path) {

    // mode == which chic
    // gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

    // loaf the fonll stuff
    TFile* fonll = new TFile((path + "/fonll.root").c_str());
    TH1F* ptHisto = (TH1F*) fonll->Get("pthisto"); 
    TH1F* etaHisto = (TH1F*) fonll->Get("etahisto"); 

    // get the graph to smear with
    TFile* smearfile = new TFile((path + "/smear12.root").c_str());
    TGraphErrors* sgraph = (TGraphErrors*)smearfile->Get("data;1");

//    //make pdfs for phi and chi mass shape
//    RooRealVar m1("m1","m1",0.6, 1.5);
//    RooRealVar m2("m2","m2",3.2, 3.6);
//    RooRealVar m3("m3","m3",3.6, 3.8);
//    RooRelBreitWigner* bw = createPhiMassPdf(m1);

//    RooRealVar lt("lt","lt",0.,1e-10);
//    RooRealVar bsGamma("bsGamma","bsGamma",1./1.51e-12);
//    RooExponential* bslifetime = new RooExponential("bslifetime","",lt,bsGamma);
//
//    // generate Bs lifetime datasets
//    RooDataSet *bslifetimedata = bslifetime->generate(RooArgSet(lt),100000) ;

//    // generate phi datasets
//    RooDataSet *phidata = bw->generate(RooArgSet(m1),100000) ;
//
//    //Double_t masses[3] = {mpi,mpi,mphi};
//    Double_t masses[3] = {0., 0.,mphi};
//    Double_t masses4[2] = {mK,mK};

//    TH1F *h1 = new TH1F("h1","h1", 100, 5150, 5550);
//    h1->Sumw2();
//    TH1F *h2 = new TH1F("mcorr","mcorr", 200, 0, 5500);
//    h2->Sumw2();
//    TH1F *dh1 = new TH1F("dh1","dh1", 100, 0., 5000.);
//    dh1->Sumw2();
//    TH1F *dh2 = new TH1F("dh2","dh2", 100, 1000., 2000.);
//    dh2->Sumw2();
//    TH1F *dh3 = new TH1F("pi","pi", 100, 2000., 4000.);
//    dh3->Sumw2();
//    TH1F *dh4 = new TH1F("phi","phi", 100, 750, 1250);
//    dh4->Sumw2();

    TRandom3 ran;

    FastDecay myDecayObject("../src/decay.txt");
    myDecayObject.setRandomGenerator(ran);
    myDecayObject.loadParentKinematics(ptHisto,etaHisto);
    myDecayObject.loadSmearGraph(sgraph);
    std::vector<int> pars;
    pars.push_back(2);
    pars.push_back(3);
    myDecayObject.addCustomParameter("qSq", FastDecay::M2, pars, true, 0., 25.);
    TH1F* arHist = new TH1F("arHist","",250,0.,25.);
    for(int i=0; i<150; ++i) {
	    arHist->SetBinContent(i+1, 1.);
    }
    myDecayObject.setAcceptRejectHist(arHist, FastDecay::M2, pars);
    myDecayObject.saveTree("myTree.root");


    int ngenerated = 0; int nselected = 0;
    for (Int_t n=0; n<nEvtToGen; ++n) {
//        double bmass;
//        if (mode == 0){
//            bmass = mBs;
//        }
//        else {
//            bmass = mBd; 
//        }
//
//        TLorentzVector Bs = genB(ran,ptHisto,etaHisto, bmass);
//        TGenPhaseSpace event;
//
//        // Generate the Bs	 
//        masses[2] = pick(phidata,ran,std::string("m1"));
//
//        if (!generateEvent(Bs,event,masses,3,ran,1000)) continue;
//        // std::cout << masses[2] << std::endl;   
//        TLorentzVector pi1 = *event.GetDecay(0); // psi
//        TLorentzVector pi2 = *event.GetDecay(1); // psi
//        TLorentzVector phi = *event.GetDecay(2); //phi
//
//        // std::cout << pi2.M() << std::endl;
//        // phi decays to kk
//        TGenPhaseSpace subevent2;  
//        if (!generateEvent(phi,subevent2,masses4,2,ran,1000)) continue;
//        TLorentzVector k1 = *subevent2.GetDecay(0);
//        TLorentzVector k2 = *subevent2.GetDecay(1);

        if (!myDecayObject.generate()) continue;
//        TLorentzVector Bs  = myDecayObject.getP(0); 
//        TLorentzVector pi1 = myDecayObject.getP(2); 
//        TLorentzVector pi2 = myDecayObject.getP(3); 
//        TLorentzVector phi = myDecayObject.getP(1); 
//        TLorentzVector k1  = myDecayObject.getP(4); 
//        TLorentzVector k2  = myDecayObject.getP(5); 
//
//        // smear the vectors
//        TLorentzVector spi1 = smearedVec(pi1,sgraph,ran);
//        TLorentzVector spi2 = smearedVec(pi2,sgraph,ran);
//        TLorentzVector sk1 =  smearedVec(k1,sgraph,ran);
//        TLorentzVector sk2 =  smearedVec(k2,sgraph,ran);
//        //           TLorentzVector sk2 = reassignMass(spi,mK);
//
//        TLorentzVector sum = spi1 + spi2 + sk1 + sk2;
//        TLorentzVector pipi = spi1 + spi2 ;
//        TLorentzVector sphi = sk1+sk2 ;
//        TLorentzVector dsum = spi1 + sk1 + sk2 ;
//        TLorentzVector dsum2 = spi2 + sk1 + sk2 ;
//
        ++ngenerated;
//        // std::cout << sum.M() << std::endl; 
//
//        //if (select(spi1,spi2,sk1,sk2,ran)) {
//        h1->Fill(sum.M()*1000);
//        dh1->Fill(dsum.M()*1000);
//        dh2->Fill(dsum2.M()*1000);
//        dh3->Fill(pipi.M()*1000);
//        dh4->Fill(sphi.M()*1000); 
        ++nselected;
//        //}
//        // std::cout << pipi.P() << std::endl;
//        
//	//get flight distance of Bs in lab frame
//        Double_t flighttime = pick(bslifetimedata,ran,std::string("lt"))*Bs.Gamma();
//        Double_t bsFD = Bs.Beta()*flighttime*3.e8;
//
//	TVector3 flight = Bs.Vect().Unit();
//	flight.SetMag(bsFD);
//
//	//smear for secondary vertex resolution (for now assume 0.02mm for x and y and 0.2mm for z)
//	//TODO should be pT dependent
//	TVector3 smearedFlight = flight + TVector3(ran.Gaus(0.,0.00002),ran.Gaus(0.,0.00002),ran.Gaus(0.,0.0002));
//	Double_t sangle = smearedFlight.Angle(sphi.Vect());
//
//	//construct the corrected mass
//	Double_t ptMiss = sphi.Rho()*TMath::Sin(sangle);
//	Double_t mPhi2  = sphi.M2();
//	Double_t mCorr = TMath::Sqrt(mPhi2 + ptMiss*ptMiss) + ptMiss;
//
//        h2->Fill(mCorr*1000);
    } //event loop	

    myDecayObject.saveHistos("myHistos.root");

//    TCanvas* can1 = new TCanvas("can","can", 800,600); 
//    h1->Fit("gaus");
//    h1->Draw("HISTO");
//    can1->Update();
//
//    TCanvas* can2 = new TCanvas("can2","can2", 800,600); 
//    dh1->Draw("HISTO");
//    can2->Update();
//
//    TCanvas* can3 = new TCanvas("can3","can3", 800,600); 
//    dh2->Draw("HISTO");
//    can3->Update();
//
//    TCanvas* can4 = new TCanvas("can4","can4", 800,600); 
//    dh3->Draw("HISTO");
//    can4->Update();
//
//    TCanvas* can5 = new TCanvas("can5","can5", 800,600); 
//    dh4->Draw("HISTO");
//    can5->Update();
//
//    TCanvas* can6 = new TCanvas("can6","can6", 800,600); 
//    h2->Draw("HISTO");
//    can6->Update();
//
    std::cout << "Generated " << ngenerated << std::endl;
    std::cout << "Selected " << nselected << std::endl;
//
//    std::stringstream outputname; outputname << "hf_signal_" << mode <<".root"; 
//
//    TFile* output = new TFile(outputname.str().c_str(),"RECREATE","output ROOT file");
//    dh1->Write(); 
//    dh2->Write(); 
//    dh3->Write(); 
//    dh4->Write(); 
//    h1->Write();  
//    h2->Write();  
//    output->Close();

}

int main(int argc, char * argv[])
{
    if (argc != 4) {
        printf("Usage: %s mode numberToGenerate pathToFiles\n", argv[0]);
        return 1;
    }
    const int mode   = atoi(argv[1]);
    const int number = atoi(argv[2]);
    const std::string path = argv[3];
    PhaseHaofei(mode, number, path);
    return 0;
}
