/*****************************************************************************
* Package: RooRarFit
 * File: $Id: RooGounarisSakurai.cc,v 1.2 2011/06/16 13:18:48 fwilson Exp $
 * Authors: A Bevan, Katherine George 
 *                                                                           
 * Gounaris-Sakurai (GS) distribution is a model of the P-wave \pi\pi         
 * scattering amplitude                                                      
 * (G.J.Gounaris and J.J.Sakurai. Phys. Rev. Lett, 21:244 (1968))             
 *
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This is an implentation of the Gounaris-Sakurai pi-pi scattering
// function for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This is an implentation of the Gounaris-Sakurai pi-pi scattering
// function for RooRarFit
// END_HTML
//

#include "Riostream.h"

#include "RooGounarisSakurai.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

//ClassImp(RooGounarisSakurai)

//---------------------------------------------------------------------------
RooGounarisSakurai::RooGounarisSakurai(const char *name, const char *title,
			 RooAbsReal& _x, RooAbsReal& _mean,
			 RooAbsReal& _width, RooAbsReal& _spin, RooAbsReal& _radius,
                         RooAbsReal& _mass_a, RooAbsReal& _mass_b) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  mean("mean","Mean",this,_mean),
  width("width","Width",this,_width),
  spin("spin","Spin",this,_spin),
  radius("radius","Form Factor Radius",this,_radius),
  mass_a("mass_a","Mass of daughter A",this,_mass_a),
  mass_b("mass_b","Mass of daughter B",this,_mass_b)
{
}

//---------------------------------------------------------------------------
RooGounarisSakurai::RooGounarisSakurai(const RooGounarisSakurai& other,
				       const char* name) : 
  RooAbsPdf(other,name), 
  x("x",this,other.x), 
  mean("mean",this,other.mean),
  width("width",this,other.width),
  spin("spin", this, other.spin), 
  radius("radius", this, other.radius), 
  mass_a("mass_a", this, other.mass_a), 
  mass_b("mass_b", this, other.mass_b) 
{
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::evaluate() const
{
  Double_t arg= x*x - mean*mean - fFunction(x);  
  Double_t gammaf = mean*Gamma();

  //  return (1 + d * width/mean)*(1 + d * width/mean) / (arg*arg + gammaf*gammaf);
  // the 1-dGamma_0/m_0 term is constant and can be ignored.  RF will deal with the
  // normalisation properly.
  if(x < (mass_a + mass_b)) return 0;
  else return x*x / (arg*arg + gammaf*gammaf);
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::Gamma() const
{
  /*
   * This is the mass (and momentum) dependent width 
   * For the GS model, there is no Blatt-Weisskopf FF ratio 
   */
  Double_t kx = KFunction((double)x);
  Double_t km = KFunction(mean);
  Double_t rk = 0.0;
  if(km!=0){
    rk = (kx/km);
    if(spin ==1){
      rk = rk*rk*rk;
    }
    if(spin ==2){
      rk = rk*rk*rk*rk*rk;
    }    
  }
  return width*(mean/x)*rk;
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::FFunction(Double_t X) const
{
  /*
   * These are the Blatt-Weisskopf form factors.  The argument 
   *  X = sqrt(s)
   */
  if(spin==0) return 1.0;
  if(spin==1) return 1.0/(1 + X*X);
  if(spin==2) return 1.0/(9 + 3*X*X + X*X*X*X); 
  return 1.0;
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::dhds() const
{
  if(mean == 0.0) return 0.0;
  Double_t k = KFunction(mean);
  if(k == 0.0) return 0.0;
  Double_t h = hFunction(mean);
  Double_t m0_2 = mean*mean;
  return h*( 1.0/(8.0*k*k) - 1.0/(2.0*m0_2)) + 1.0/(2*M_PI*m0_2);
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::hFunction(Double_t X) const
{
  /*
   * Argument X = sqrt(s)
   */
  if(X == 0.0) return 0.0;
  // assume that the pion mass is the average of the two pion masses used.
  Double_t mpi  = 0.5*(mass_a + mass_b); 
  Double_t k = KFunction(X);

  Double_t theLog =  log( (X + 2.0*k) / (2.0* mpi) );
  return (2.0*k/(M_PI*X))*theLog;
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::fFunction(Double_t X) const
{
  /*
   * Argument is X = sqrt(s)
   */
  Double_t grad = dhds();
  Double_t h_s  = hFunction(X);
  Double_t h_m0 = hFunction(mean);
  Double_t k_s  = KFunction(X);
  Double_t k_m0 = KFunction(mean);
  if(k_m0 == 0.0) return 0.0;
  Double_t k2_m0 = k_m0*k_m0;
  Double_t k3_m0 = k2_m0*k_m0;
  Double_t mean2 = mean*mean;

  Double_t func = k_s*k_s*(h_s - h_m0) + (mean2 - X*X)*k2_m0*grad;
  return func * width*mean2 / k3_m0;
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::KFunction(Double_t X) const
{
  /*
   * This is momentum calculation
   * Argument X = sqrt(s)
   */
  if(X==0) return 0;
  Double_t massone = sqrt(1 - (mass_a + mass_b)*(mass_a + mass_b)/(X*X));
  Double_t masstwo = sqrt(1 - (mass_a - mass_b)*(mass_a - mass_b)/(X*X));
  return X/2.0 * massone * masstwo;
}

//---------------------------------------------------------------------------
Double_t RooGounarisSakurai::dFunction() const
{
  // assume that the pion mass is the average of the two pion masses used.
  Double_t mpi  = 0.5*(mass_a + mass_b); 
  Double_t mpi2 = mpi*mpi; 
  Double_t kfunc = KFunction(mean); 
  Double_t kfunc2 = kfunc*kfunc;

  Double_t logCoeff = (3*mpi2)/(M_PI*kfunc2);
  Double_t one = (mean + 2*kfunc)/(2*mpi);
  Double_t two = mean/(2*M_PI*kfunc);
  Double_t three = (mpi2*mean)/(M_PI*kfunc2*kfunc);

  return (logCoeff * log(one) + two - three);
}




