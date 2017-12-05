#ifndef ROO_RELBREITWIGNER
#define ROO_RELBREITWIGNER

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooRelBreitWigner : public RooAbsPdf {

public:

  inline RooRelBreitWigner() {};

  RooRelBreitWigner(const char *name, 
		    const char *title,
		    RooAbsReal& _x, 
		    RooAbsReal& _mean, 
		    RooAbsReal& _width,
		    RooAbsReal& _radius,
		    RooAbsReal& _mass_a,
		    RooAbsReal& _mass_b,
		    RooAbsReal& _spin);
  RooRelBreitWigner(const RooRelBreitWigner& other, 
		    const char* name=0) ;
  virtual TObject* clone(const char* newname) const { 
    return new RooRelBreitWigner(*this,newname); 
  }
  inline virtual ~RooRelBreitWigner() { }

protected:

  RooRealProxy x ;      // observable (mass)
  RooRealProxy mean ;    // mean
  RooRealProxy width ;  // width
  RooRealProxy radius ; // meson radius parameter
  RooRealProxy mass_a ; // mass of first daughter from decay
  RooRealProxy mass_b ; // mass of second daughter from decay
  RooRealProxy spin ; // spin (= 0,1,2)
  
  Double_t evaluate() const ;

private:
  Double_t getWidth() const;
  Double_t getQterm(Double_t q, Double_t q0) const;
  Double_t getFF(Double_t q) const;
  Double_t getQ(Double_t mass) const;

  //ClassDef(RooRelBreitWigner,0) // Relativistic Breit Wigner PDF
};

#endif
