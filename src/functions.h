#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <string>

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "RooRelBreitWigner.h"
#include "RooGounarisSakurai.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

const double me      = 0.000510998928;
const double mmu     = 0.1134289267;
const double mtau    = 1.77686;
const double mpi0    = 0.1349766;
const double mpi     = 0.13957;
const double meta    = 0.547862;
const double mrho0   = 0.77526;
const double mrho    = 0.77511;
const double momega  = 0.78265;
const double mK0     = 0.497611;
const double mK      = 0.493677;
const double mKstar  = 0.89166;
const double mKstar0 = 0.89581;
const double metapr  = 0.95778;
const double mphi    = 1.0197;
const double mDplus  = 1.86961;
const double mD0     = 1.86484;
const double mDs     = 1.9683;
const double mDstar  = 2.01027;
const double mD0star = 2.00697;
const double mDsstar = 2.1121;
const double metac   = 2.9836;
const double mchic0  = 3.41475;
const double mJpsi   = 3.096916;
const double mchic   = 3.51066;
const double mPsi    = 3.686093;
const double mchic2  = 3.55620;
const double mBd     = 5.27958;
const double mB      = 5.27929;
const double mBs     = 5.36677;
const double mBc     = 6.2751;
const double mB0star = 5.32483;
const double mBstar  = 5.32483;
const double mBsstar = 5.4154;
const double metab   = 9.3980;
const double mUp1S   = 9.46030;
const double mUp2S   = 10.02326;
const double mUp3S   = 10.3552;
const double mUp4S   = 10.5794;

const double mn      = 0.939565379;
const double mp      = 0.938272046;
const double mL      = 1.115683;
const double mSm     = 1.197449;
const double mS0     = 1.192642;
const double mSp     = 1.18937;
const double mLc     = 2.28646;
const double mSc0    = 2.45375;
const double mScp    = 2.4529;
const double mScpp   = 2.45397;
const double mLb     = 5.61951;
const double mSbm    = 5.8155;
const double mSbp    = 5.8113;

bool generateEvent(TLorentzVector& head, TGenPhaseSpace& event, double* masses , int np, TRandom& ran, int m_maxgen );

double pick(RooDataSet * data ,TRandom& ran, std::string var_name);

const double barrierFactor = 3;

std::string varname(std::string header, std::string var);

RooGounarisSakurai* rooGS(RooRealVar& m, double mean = 0.77511 , double gamma = 0.1491,
        double thespin = 1, double m1 = 0.13957,double m2 = 0.1349766 , std::string name = "rhoplus"); 


RooRelBreitWigner* rooBW(RooRealVar& m,  double mean = 1.019461 , double gamma = 0.001491,
        double thespin = 1, double m1 =0.493677,double m2 = 0.493677 , std::string name = "phi");

RooGounarisSakurai* createRhoPlus(RooRealVar& m, std::string name = "rhoplus");

RooRelBreitWigner* createPhiMassPdf(RooRealVar& m, std::string name = "phi");

RooRelBreitWigner* createKstarMassPdf(RooRealVar& m, std::string name = "kstar");

RooRelBreitWigner* createChi0MassPdf(RooRealVar& m, std::string name = "chic0");

RooRelBreitWigner* createChi1MassPdf(RooRealVar& m, std::string name = "chic1");

RooRelBreitWigner* createChi2MassPdf(RooRealVar& m, std::string name = "chic2");

RooRelBreitWigner* createpsi2MassPdf(RooRealVar& m, std::string name = "chic2");

TLorentzVector genB(TRandom ran, TH1F* ptHisto, TH1F* etaHisto, double m);

double resSlope(double p);

TLorentzVector toFourVector(const TVector3& vec, double m);

TLorentzVector reassignMass(const TLorentzVector& vec, double mass);

TVector3 smearedVec(double plus_x, double plus_y, double plus_z, TGraphErrors* dGraph, TRandom& ran);

TLorentzVector smearedVec(TLorentzVector& vec, TGraphErrors* dGraph, TRandom& ran);

bool inAcceptance(TLorentzVector& vec);

const double zC = 5.4;
const double ptkick= 1.2;
const double zTracker = 9.5; // end of tracker...
const double xSizeTracker = 9.5*0.3;
const double ySizeTracker = 9.5*0.25;
const double xMinTracker = 9.5*0.001;
const double yMinTracker = 9.5*0.001;


TLorentzVector magnetKick(TLorentzVector& vec, double charge);

bool inDownstream(TLorentzVector& vec, int charge);

int pdgCode(TString part);

double getMass(int pdgCode);

#endif

/// common functions
