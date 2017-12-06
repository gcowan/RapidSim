// Written by Philip Ilten and Mike Williams, 20/06/2017.
// Toolset to check if a secondary vertex is consistent with having
// originated from the VELO material. See README for an overview.
#ifndef VELOMATERIAL_H
#define VELOMATERIAL_H
#include <limits>
#include <algorithm>
#include "TFile.h"
#include "TKey.h"
#include "TObjString.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TError.h"
#include "TGraph.h"
#include "TSpline.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
using namespace std;

// Pre-declare material class for sensor configuration.
class FoilMaterial;

//=============================================================================
/*
  Class to access the parameterized VELO module material.
*/
class ModuleMaterial {
  
public:

  //===========================================================================
  /*
    Class to hold the module sensor geometry.

    Details on the module sensor geometry description can be found in
    the DBASE source:
  
    DBASE/Det/XmlDDDB/DDDB/Velo/Geom/LogVol/GenericModule.xml
    DBASE/Det/XmlDDDB/DDDB/Velo/GeomParam/VeloSensorParam.xml
    
    Here, the module geometry was taken from the EDMS design
    documents with number 401568 v.4 that can be found at the
    following link.
    
    https://edms.cern.ch/ui/#!master/navigator/document?D:1100625139:
    1100625139:subDocs
  */
  class Sensor {
    
  public:
    
    // Configurable Members.
    TF1 fncl, fncu; ///< Functions defining the lower/upper sensor edges.
    double z;       ///< Sensor z-position.
    int half;       ///< Lower/upper VELO half.

    // Methods.
    /// Default constructor.
    Sensor();
    
    /*
      Primary constructor.
      \param ifile: TFile containing all the needed parameters.
      \param mod:   string of the module to build, e.g. "module_01".
      \param sns:   string of the sensor to build from the module, e.g. "0".
    */
    Sensor(TFile &file, string mod, string sns);

    /// Operator for sorting by the given z-value.
    bool operator<(const Sensor& sns) const;

    /// Configure the tool.
    void config(double xc, double yc, bool rel = true);

    /// Configure the tool.
    void config(const TF1 &fnc);

    /*
      Return if an x,y-point falls within a module sensor.
      \param [x,y]: x,y-position to evaluate.
    */
    bool inside(double x, double y);

    /*
      Return the x-position for the edge of the sensor.
      \param y:    y-position to evaluate.
      \param edge: -/+1 for the lower/upper edge, 0 for the inside edge.
    */
    double x(double y, int edge = 0);

    /*
      Return a parameter for the sensor function.
      \param par:  parameter index.
      \param edge: -/+1 for the lower/upper edge, 0 for the inside edge.
    */
    double p(unsigned int par, int edge = 0);
  };

  //===========================================================================
  /*
    Class to calculate the probability integral for a vertex to be
    produced from VELO modules.

    The probability density function for the vertex is assumed to be a
    3D Gaussian and the integration is performed in the x,y-plane as
    the z-direction pre-factor remains constant.
  */
  class Integral {
    
    // Members.
    vector<Sensor> *mzs;                 ///< Module sensors.
    Sensor         *geo;                 ///< Sensor geometry for integration.
    TF2 fnc;                             ///< Function used for integration.
    ROOT::Math::WrappedMultiTF1 wrp;     ///< Wrapper for the function.
    ROOT::Math::IntegratorMultiDim *inr; ///< Integrator for the function.
    double vx, vy, vz, dvx, dvy, dvz;    ///< Vertex position and uncertainty.

    /*
      Return the probability integral for a given half.
      \param halfIn: -/+1 for the lower/upper VELO half in x.
    */
    double integral(int half);

  public:  

    // Configurable members.
    double dvmin; ///< Minimum allowed vertex uncertainty.
    double sigma; ///< Number of standard deviations to calculate over.
    double width; ///< Width of the module in mm.
    double tol;   ///< Integration tolerance.
    /// Integration algorithm: adaptive, VEGAS, MISER, or plain.
    ROOT::Math::IntegrationMultiDim::Type alg; 

    // Methods.
    /*
      Default constructor.
      \param mzsIn: module sensors.
    */
    Integral(vector<Sensor> *mzs);

    ///< Destructor.
    ~Integral(); 

    /*
      Return the probability function for a given half. The passed
      coordinates are not in mm but are in sigma from the origin of a
      2D normal distribution.
      \param t: the variables sx and sy passed as an array of length 2.
    */
    double operator()(double *t, double *);
    
    /// Configure the tool.
    void config(double dvminIn, double sigmaIn, double mwIn,
		double tolIn, ROOT::Math::IntegrationMultiDim::Type algIn);

    /// Return the probability integral.
    double integral(double vxIn, double vyIn, double vzIn, double dvxIn,
		    double dvyIn, double dvzIn, int halfIn);
  };

  //===========================================================================
  /*
    Class to calculate the minimum distance in uncertainty space
    between a vertex and the VELO modules.

    First the lower and upper modules which minimize the distance in z
    are found. This is the same in normal or uncertainty space. Next,
    the minimum transverse distance is found for each module in
    uncertainty space. The minimum distance between the the lower and
    upper modules is returned.
  */
  class Distance {
    
    // Members.
    vector<Sensor> *mzs;              ///< Module sensors.
    Sensor         *geo;              ///< Sensor geometry for intersection.
    TF1 fnc;                          ///< Wrapped function to minimize.
    double vx, vy, vz, dvx, dvy, dvz; ///< Vertex position and uncertianty.
  
    /*
     Return the distance for a given half.
     \param halfIn: -/+1 for the lower/upper VELO half in x.
     \param mzs:    module z-values.
     \param s:      distance vector for updating.
    */
    void distance(int half, TLorentzVector &s);
 
  public:
    
    // Configurable members.
    double dvmin; ///< Minimum allowed vertex uncertainty. 
    double tol;   ///< Minimization tolerance.
    int    itr;   ///< Minimization iterations.
    bool   logt;  ///< Flag to use a log scale in minimization.

    // Methods.
    /*
      Default constructor.
      \param mzsIn: module sensors.
    */ 
   Distance(vector<Sensor> *mzsIn);
    
    /*
      Return the transverse distance.
      \param t: the variable sy passed as an array of length 1.
    */
    double operator()(double *t, double *);

    /// Configure the tool.
    void config(double dvminIn, double tolIn, int itrIn, bool logtIn);

    /// Return the distance.
    TLorentzVector distance(double vxIn, double vyIn, double vzIn,
			    double dvxIn, double dvyIn, double dvzIn,
			    int halfIn);
  };

  //===========================================================================
private:

  // Members.
  vector<Sensor> mzs; ///< Module sensors.
  Integral  ingr;     ///< Internal integral tool.
  Distance  dist;     ///< Internal distance tool.

public:

  // Methods.
  /*
    Default constructor.
    \param pars: name of the fit parameters file.
  */
  ModuleMaterial(string pars = "dat/pars.root");

  /*
    Return the index for the sensor nearest a given z-position.
    \param z:    z-position to evaluate.
    \param half: -/+1 for the lower/upper VELO half in x, 0 for both.
   */
  int sensor(double z, int half = 0);

  /*
    Return if an x,y-point falls within a module sensor.
    \param idx: index of the module sensor.
    \param x:   x-position to evaluate.
    \param y:   y-position to evaluate.
  */
  bool inside(unsigned int idx, double x, double y);

  /*
    Return the x-position for the edge of a module sensor.
    \param idx:  index of the module sensor.
    \param y:    y-position to evaluate.
    \param edge: -/+1 for the lower/upper edge, 0 for the inside edge.
   */
  double x(unsigned int idx, double y, int edge = 0);

  /*
    Return the z-position for a module sensor.
    \param idx:  index of the module sensor.
   */
  double z(unsigned int idx);

  /*
    Return a parameter for a module sensor function.
    \param idx:  index of the module sensor.
    \param par:  parameter index.
    \param edge: -/+1 for the lower/upper edge, 0 for the inside edge.
  */
  double p(unsigned int idx, int par, int edge = 0);

  /*
    Configure the technical sensor parameters. This method either sets
    the absolute x,y-position for the center of each sensor or shifts
    the center from the nominal position. This can be applied to
    either or both halves of the VELO sensors.
    \param xc:   x-position of the sensor center in mm.
    \param yc:   y-position of the sensor center in mm.
    \param rel:  if true, shift the sensors from their nominal positions.
    \param half: -/+1 for the lower/upper VELO half in x, 0 for both.
   */ 
  void configSensor(double xc = 0, double yc = 0, bool rel = true, 
		    int half = 0);

   /*
    Configure the technical sensor parameters. The sensor edge is set
    as the foil edge at the given z-position for each sensor.
    \param foil: the foil material map.
   */ 
  void configSensor(FoilMaterial &foil);

  /*
    Configure the technical probability integral parameters.
    \param dvmin: minimum allowed vertex uncertainty.
    \param sigma: number of standard deviations to calculate over.
    \param width: width of the module in mm.
    \param tol:   integration tolerance.
    \param alg:   integration algorithm: adaptive, VEGAS, MISER, or plain.
  */
  void configIntegral(double dvmin = 0.25, double sigma = 5,
		      double width = 1, double tol = 0.01,
		      ROOT::Math::IntegrationMultiDim::Type alg = 
		      ROOT::Math::IntegrationMultiDim::kADAPTIVE);

  /*
    Configure the technical distance parameters.
    \param dvmin: minimum allowed vertex uncertainty.
    \param tol:   minimization tolerance.		  
    \param itr:   minimization iterations		  
    \param logt:  flag to use a log scale in minimization.
  */
  void configDistance(double dvmin = 0.25, double tol = 1e-10, int itr = 100,
		      bool logt = false);

  /*
    Return the probability integral.
    \param v[x,y,z]:    x,y,z-positions of the vertex.
    \param d[vx,vy,vz]: uncertainties on the x,y,z-positions of the vertex.
    \param half:        -/+1 for the lower/upper VELO half in x, 0 for both.
  */
  double integral(double vx, double vy, double vz,
		  double dvx, double dvy, double dvz, int half);

  /*
    Return the intersection. A 4-vector of the intersection point is
    returned, where the T (or E) component gives the VELO half (-/+1
    for lower/upper) or 0 if no intersection is found.
    \param v[x,y,z]: x,y,z-positions of the vertex.
    \param f[x,y,z]: flight direction.
    \param half:     -/+1 for the lower/upper VELO half in x, 0 for both.
  */
  TLorentzVector intersect(double vx, double vy, double vz,
			   double fx, double fy, double fz, int half = 0);

  /*
    Return the distance. The point of minimum distance on the material
    in uncertainty space is returned where the T (or E) component
    gives the distance.
    \param v[x,y,z]:    x,y,z-positions of the vertex.
    \param d[vx,vy,vz]: x,y,z-uncertainties of the vertex position.
    \param half:        -/+1 for the lower/upper VELO half in x, 0 for both.
  */
  TLorentzVector distance(double vx, double vy, double vz,
			  double dvx = 1, double dvy = 1, double dvz = 1,
			  int half = 0);
};

//=============================================================================
/*
  Class to access the parameterized VELO foil material.
*/
class FoilMaterial {
  friend class ModuleMaterial;
  
public:

  //===========================================================================
  /*
    Class to hold the foil segment geometry.

    This class is built from the parameters and functions provided in
    the input parameter file passed to the primary constructor. The
    class holds the function for the x,y-slice as well as the
    parameters for this slice as a function of z-position.
  */
  class Segment {
    friend class ModuleMaterial;

    // Members.
    vector<TF1>      fncs; ///< Vector of functions for the slice parameters.
    vector<TSpline3> spls; ///< Vector of splines for the slice parameters.
    vector<int>      mths; ///< Vector of methods for parameter evaluation.
    TF1 fnc;               ///< The slice function.
    double zmin, zmax;     ///< Range of validity in z-position for the segment.
    
  public:
    
    // Methods.
    /// Default constructor.
    Segment();
    
    /*
      Primary constructor.
      \param ifile: TFile containing all the needed parameters.
      \param pre:   string of the foil segment to build, e.g. "foil_l_c".
    */
    Segment(TFile &file, string pre);
    
    /// Configure the tool.
    void config(int par, int mth);

    /*
      Return if a z-position is in the segment.
      \param z: z-position to evaluate.
    */
    bool valid(double &z);

    /*
      Return the x-position for the segment.
      \param [y,z]: y,z-position to evaluate.
    */
    double x(double y, double z);

    /*
      Return a parameter for the segment function.
      \param par: parameter index.
      \param z:   z-position to evaluate.
    */
    double p(unsigned int par, double z);
  };
  
  //===========================================================================
  /*
    Class to calculate the probability integral for a vertex to be
    produced from the VELO foil.

    The probability density function for the vertex is assumed to be a
    3D Gaussian and the integration is performed in the y,z-plane as
    the x-direction pre-factor remains constant.
  */
  class Integral {
    
    // Members.
    FoilMaterial *foil;                  ///< Foil geometry.
    TF2 fnc;                             ///< Function used for integration.
    ROOT::Math::WrappedMultiTF1 wrp;     ///< Wrapper for the function.
    ROOT::Math::IntegratorMultiDim *inr; ///< Integrator for the function.
    double vx, vy, vz, dvx, dvy, dvz;    ///< Vertex position and uncertainty.
    int half;                            ///< Lower/upper VELO half to evaluate.

    /*
      Return the probability integral in one dimension.
      \param v:  vertex position.
      \param dv: vertex position uncertainty.
      \param t:  material position.     
      \param dt: material position width.
    */
    double integral(double v, double dv, double t, double dt);
    
  public:  

    // Configurable members.
    double dvmin; ///< Minimum allowed vertex uncertainty.
    double sigma; ///< Number of standard deviations to calculate over.
    double step;  ///< Step distance in sigma.
    double width; ///< Width of the foil in mm.
    double tol;   ///< Integration tolerance.
    /// Integration algorithm: adaptive, VEGAS, MISER, or plain.
    ROOT::Math::IntegrationMultiDim::Type alg; 

    // Methods.
    /*
      Default constructor.
      \param foilIn: foil geometry.
    */
    Integral(FoilMaterial *foilIn);

    ///< Destructor.
    ~Integral(); 

    /*
      Return the probability function for a given half. The passed
      coordinates are not in mm but are in sigma from the origin of a
      2D normal distribution.
      \param t: the variables sy and sz passed as an array of length 2.
    */
    double operator()(double *t, double *);

    /// Configure the tool.
    void config(double dvminIn, double sigmaIn, double stepIn, double widthIn,
		double tolIn, ROOT::Math::IntegrationMultiDim::Type algIn);

    /// Return the probability integral.
    double integral(double vxIn, double vyIn, double vzIn, double dvxIn,
		    double dvyIn, double dvzIn, int methodIn);
  };

  //===========================================================================
  /*
    Class to calculate the intersection of a particle, defined with a
    vertex and flight direction, with the VELO foil.
  */
  class Intersect {
    
    // Members.
    FoilMaterial *foil;               ///< Foil geometry.
    TF1 fnc;                          ///< Wrapped function to minimize.
    double vx, vy, vz, fx, fy, fz;    ///< Vertex position and flight direction.
    int half;                         ///< Lower/upper VELO half to evalulate.

    /*
      Return the intersection for a given half and range.
      \param halfIn: -/+1 for the lower/upper VELO half in x.
      \param t:      point along flight path to evaluate.
    */
    double intersect(int halfIn, double t);

    /*
      Check if a point is within valid bounds.
      \param t: distance along flight path.
    */
    bool check(double t);

  public:
    
    // Configurable members.
    double step; ///< Step distance in mm.
    double tol;  ///< Minimization tolerance.
    int    itr;  ///< Minimization iterations.
    bool   logt; ///< Flag to use a log scale in minimization.

    // Methods.
    /*
      Default constructor.
      \param foilIn: foil geometry.
    */
    Intersect(FoilMaterial *foilIn);
    
    /*
      Return the distance from the foil along the flight direction.
      \param t: the variable t passed as an array of length 1.
    */
    double operator()(double *t, double *);

    /// Configure the tool.
    void config(double stepIn, double tolIn, int itrIn, bool logtIn);

    /// Return the intersection.
    TLorentzVector intersect(double vxIn, double vyIn, double vzIn,
			     double fxIn, double fyIn, double fzIn,
			     int halfIn);
  };

  //===========================================================================
  /// Class to calculate the distance between a vertex and the foil.
  class Distance {
    
    // Members.
    FoilMaterial *foil;               ///< Representation of the foil.
    ROOT::Math::Functor fnc;          ///< Wrapped function to minimize.
    ROOT::Math::Minimizer *min;       ///< Minimizer.
    double vx, vy, vz, dvx, dvy, dvz; ///< Vertex position and uncertainty.
    int half;                         ///< Lower/upper VELO half to evaluate.
    int level;                        ///< Saved ROOT reporting level.

    /*
      Return the distance for a given half.
      \param halfIn: -1 for the lower half and +1 for the upper foil half.
      \param method: distance method.
      \param dir:    minimizing z-direction: 0 both, -/+1 backwards/forward.
      \param s:      vector containt
      \param n:      number of all distances.
      \param u:      sum of all distances.
    */
    void distance(int halfIn, int method, int dir, TLorentzVector &s,
		  double &n, double &u);
    
  public:
    
    // Configurable members.
    double dvmin; ///< Minimum allowed vertex uncertainty.
    double dvmax; ///< Maximum vertex uncertainty for minimization limits.
    double sigma; ///< Number of standard deviations to calculate over.
    double step;  ///< Minimization step size along z-direction in mm.
    double tol;   ///< Minimization tolerance.
    double itr;   ///< Minimization iterations.
    string lib;   ///< Minimization library, see alg for options.
    /*
      Minimization algorithm. The available algorithms for each library are:
      \param Minuit2:     Migrad, Simplex, Combined, Scan, Fumili
      \param GSLMultiMin: ConjugateFR, ConjugatePR, BFGS, BFGS2, 
                          SteepestDescent
      \param GSLMultiFit: 
      \param GSLSimAn: 
    */
    string alg;
    
    // Methods.
    /*
      Default constructor.
      \param foilIn: foil geometry.
    */
    Distance(FoilMaterial *foilIn);
    
    /*
      Return the distance from the foil in uncertainty space.
      \param t: the variables y and z passed as an array of length 2.
    */
    double operator()(const double *t);

    /// Configure the tool.
    void config(double dvminIn, double dvmax, double sigmaIn, double stepIn,
		double tolIn, int itrIn, string libIn, string algIn);

    /// Return the distance.
    TLorentzVector distance(double vxIn, double vyIn, double vzIn,
			    double dvxIn, double dvyIn, double dvzIn,
			    int halfIn, int method, int dir);
  };

  //===========================================================================
protected:

  // Members.
  vector<Segment> flzs, fuzs; ///< Lower/upper foil segment functions.
  Integral  ingr;             ///< Internal integral tool.
  Intersect insc;             ///< Internal intersection tool.
  Distance  dist;             ///< Internal distance tool.

public:

  // Methods.
  /*
    Default constructor.
    \param pars: name of the fit parameters file.
  */
  FoilMaterial(string pars = "dat/pars.root");
  
  /*
    Return the x-position of the foil.
    \param [y,z]: y,z-position to evaluate.
    \param half:  -/+1 for the lower/upper VELO half in x.
  */
  double x(double y, double z, int half);

  /*
    Return a parameter for the foil function.
    \param par:  parameter index.
    \param z:    z-position to evaluate.
    \param half: -/+1 for the lower/upper VELO half in x.
  */
  double p(unsigned int par, double z, int half);

  /*
    Configure the technical segment parameters.
    \param half: -/+1 for the lower/upper VELO half in x, 0 for both.
    \param seg:  name of the segment: "b", "c", "t", "f". If "a" all segments 
                 are configured.
    \param par:  parameter to configure, -1 for all.
    \param mth:  method for parameter evaluation: 0 for fitted function and
                 1 for interpolation.
  */
  void configSegment(int half = 0, string seg = "a", int par = -1, int mth = 0);

  /*
    Configure the technical probability integral parameters.
    \param dvmin: minimum allowed vertex uncertainty.
    \param sigma: number of standard deviations to calculate over.
    \param step:  integration step size along z-direction in sigma. 
    \param width: width of the foil in mm.
    \param tol:   integration tolerance.
    \param alg:   integration algorithm.
  */
  void configIntegral(double dvmin = 0.25, double sigma = 5, double step = 1, 
		      double width = 1, double tol = 1e-2,
		      ROOT::Math::IntegrationMultiDim::Type alg = 
		      ROOT::Math::IntegrationMultiDim::kADAPTIVE);

  /*
    Configure the technical intersection parameters.
    \param step: step distance in mm.
    \param tol:  minimization tolerance.
    \param itr:  minimization iterations
    \param logt: flag to use a log scale in minimization.
  */
  void configIntersect(double step = 1, double tol = 1e-10,
		       int itr = 100, bool logt = false);

  /*
    Configure the technical distance parameters.
    \param dvmin: minimum allowed vertex uncertainty.
    \param dvmax: maximum vertex uncertainty for minimization limits.
    \param sigma: number of standard deviations to calculate over.
    \param step:  minimization step size along z-direction in mm.
    \param tol:   minimization tolerance.
    \param itr:   minimization iterations.
    \param lib:   library used for minimization.
    \param alg:   algorithm used for minimization.
  */
  void configDistance(double dvmin = 0.25, double dvmax = 100,
		      double sigma = 5, double step = 1,
		      double tol =  1e-2, int itr = 100, 
		      string lib = "Minuit2", string alg = "Migrad");

  /*
    Return the probability integral.
    \param v[x,y,z]:    x,y,z-positions of the vertex.
    \param d[vx,vy,vz]: uncertainties on the x,y,z-positions of the vertex.
  */
  double integral(double vx, double vy, double vz,
		  double dvx, double dvy, double dvz, int method = 0);

  /*
    Return the intersection. A 4-vector of the intersection point
    is returned, where the T (or E) component gives the VELO half
    (-/+1 for lower/upper) or 0 if no intersection is found.
    \param v[x,y,z]: x,y,z-positions of the vertex.
    \param f[x,y,z]: flight direction.
    \param half:     -/+1 for the lower/upper VELO half in x, 0 for both.
  */
  TLorentzVector intersect(double vx, double vy, double vz, double fx,
			   double fy, double fz, int half = 0);

  /*
    Return the distance. The point of minimum distance on the material
    in uncertainty space is returned where the T (or E) component
    gives the distance.
    \param v[x,y,z]:    x,y,z-positions of the vertex.
    \param d[vx,vy,vz]: x,y,z-uncertainties of the vertex position.
    \param half:        -/+1 for the lower/upper VELO half in x, 0 for both.
    \param method:      distance method.
    \param dir:         minimizing z-direction: 0 both, -/+1 backwards/forward.
    Three distance methods are implemented. (1) The minimum
    distance. (2) The inverse of the sum of lower and upper minimum
    distance reciprocals. (3) Same as the previous method, but now
    rather than the minimum distances, it is summed over all distances.
  */
  TLorentzVector distance(double vx, double vy, double vz,
			  double dvx = 1, double dvy = 1, double dvz = 1,
			  int half = 0, int method = 0, int dir = 0);
};

//=============================================================================
/*
  Class to access the parameterized VELO material, both module and
  foil, as implemented in the initial inclusive dimuon dark photon
  analysis, LHCb-ANA-2017-027.
*/
class VeloMaterial {

public:

  /*
    Default constructor.
    \param pars: name of the fit parameters file.
  */
  VeloMaterial(string pars = "dat/pars.root");

  /*
    Return the distance. This is the harmonic mean of the following
    six distances: upper module, lower module, forward upper foil,
    forward lower foil, backward upper foil, backward lower foil.
    \param v[x,y,z]:    x,y,z-positions of the vertex.
    \param d[vx,vy,vz]: x,y,z-uncertainties of the vertex position.
  */
  double distance(double vx, double vy, double vz,
		  double dvx, double dvy, double dvz);
  
  /*
    Return true if the flight direction of a vertex intersects a foil
    tip, i.e. intersects the foil within 1 mm of where the modules are
    located.
    \param v[x,y,z]: x,y,z-positions of the vertex.
    \param f[x,y,z]: flight direction.
   */
  bool tip(double vx, double vy, double vz, double fx, double fy, double fz);
  
  /*
    Return true if the expected first hit is missed, given a flight direction
    and vertex. The active area of the sensors begins at 8.170 mm.
    \param v[x,y,z]: x,y,z-positions of the vertex.
    \param f[x,y,z]: flight direction.
    \param mz: z-position of the module containing the first hit.
  */
  bool miss(double vx, double vy, double vz, double fx, double fy, double fz,
	    double mz);
  
 private:

  // Members.
  FoilMaterial   foil; ///< Foil material tool.
  ModuleMaterial modf; ///< Module material tool extended to the foil.
  ModuleMaterial modm; ///< Module material tool with standard geometry.
  
};

#endif // VELOMATERIAL_H
