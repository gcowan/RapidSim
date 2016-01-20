/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 * File: $Id: RooGounarisSakurai.rdl,v 1.2 2011/06/16 13:18:48 fwilson Exp $
 * Authors: Katherine George - University of Liverpool                       *
 * 
 *****************************************************************************/
#ifndef ROO_GOUNARISSAKURAI
#define ROO_GOUNARISSAKURAI

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooGounarisSakurai : public RooAbsPdf {

public:
  RooGounarisSakurai(const char *name, 
		     const char *title,
		     RooAbsReal& _x, 
		     RooAbsReal& _mean, 
		     RooAbsReal& _width, 
		     RooAbsReal& _spin, 
		     RooAbsReal& _radius,
		     RooAbsReal& _mass_a, 
		     RooAbsReal& _mass_b);
  RooGounarisSakurai(const RooGounarisSakurai& other, 
		     const char* name=0) ;
  virtual TObject* clone(const char* newname) const { 
    return new RooGounarisSakurai(*this,newname); 
  }
  inline virtual ~RooGounarisSakurai() { }

 protected:

  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy spin;
  RooRealProxy radius;
  RooRealProxy mass_a;
  RooRealProxy mass_b;
 
  Double_t evaluate() const ;

private:

  Double_t Gamma() const;
  Double_t KFunction(Double_t X) const;
  Double_t FFunction(Double_t X) const;
  Double_t dFunction() const;
  Double_t hFunction(Double_t X) const;
  Double_t dhds() const;
  Double_t fFunction(Double_t X) const;

  //ClassDef(RooGounarisSakurai,0) // Gounaris Sakurai PDF
};

#endif
