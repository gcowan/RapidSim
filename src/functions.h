#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "TLorentzVector.h"
#include "TMath.h"

TLorentzVector toFourVector(const TVector3& vec, double m);

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

#endif

/// common functions
