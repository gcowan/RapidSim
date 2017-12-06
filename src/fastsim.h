#ifndef FASTSIM_H
#define FASTSIM_H
#include "velo.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip>
#define PI 3.14159265358979323846

//=============================================================================
/*
  Class to hold a particle. TO-DO: expand to hold all needed info.
*/
class Particle : public TLorentzVector {
public:
  /*
    Constructor.
  */
  Particle(const double &px, const double &py, const double &pz,
	   const double &e, int idIn = 0);
  
  // Members.
  TLorentzVector x; ///< Most recent interaction point of the particle.
  int id;           ///< Particle ID.
  int nvelo;        ///< Number of VELO hits.

  // Friends.
  friend ostream &operator<<(ostream &s, const Particle &p);
};

//=============================================================================
/*
  Class to hold an event. TO-DO: expand to hold all needed info.
*/
class Event : public vector<Particle> {
public:
  // Friends.
  friend ostream &operator<<(ostream&, const Event&);
};

//=============================================================================
/*
  Class to provide constants for a material. All distance units are in
  mm. Energy units are not specified.
*/
class Material {
public:

  /*
    Primary constructor.
    \param name: symbol of the material element, e.g. "Si" or "Al".
  */
  Material(string name);

  // Constant members.
  const double alpha = 0.00729735; ///< Fine structure constant.
  const double re    = 2.8179e-12; ///< Classical radius of the electron.

  // Members.
  double z;  ///< Atomic number.
  double na; ///< Number of atoms per volume.
  double fe; ///< Elastic form factor.
  double fi; ///< Inelastic form factor.
  double fc; ///< Form factor Coulomb correction.
};

//=============================================================================
/*
  Base class for physics processes.
*/
class Process {
public:
  /*
    Base destructor.
  */
  virtual ~Process() {;}

  /*
    Perform a physics process on a given particle in the event record.
    \param e: the event record.
    \param i: index of the particle to process.
    \param m: material inducing the process.
    \param w: width of the material in mm.
  */
  virtual void process(Event &e, const int &i, const Material &m, 
		       const double &w) = 0;
};

//=============================================================================
/*
  Class to perform Bremsstralung radiation for a charged particle
  passing through material.
*/
class ProcessBrem : public Process {
public:

  /*
    Constructor.
    \param rndIn: random number generator.
  */
  ProcessBrem();

  /*
    Add bremsstrahlung to the event record for a given particle. This
    method calculates the probability for a non-emission, determines
    the distance for this, and then generates a photon at that
    distance.
  */
  void process(Event &e, const int &i, const Material &m, const double &w);

  /*
    Sample an emitted photon energy from the differential cross-section.
    \param m:  material inducing the bremsstrahlung.
    \param ec: energy of the charged radiator.
  */
  double energy(const Material &m, const double &ec);

  /*
    Sample the polar angle theta of the emitted photon, with respect
    to the charged radiator.
    \param mc: mass of the charged radiator.
    \param ec: energy of the charged radiator.
  */
  double theta(const double &mc, const double &ec);

  /*
    Sample the azimuthal angle phi of the emitted photon, with respect
    to the charged radiator. This is just a uniform distribution
    between 0 and 2pi.
  */
  double phi();

  /*
    Calculate the differential cross-section in emitted photon energy.
    \param m:  material inducing the bremsstrahlung.
    \param eg: energy of the emitted photon.
    \param ec: energy of the charged radiator.
  */
  double dsigma(const Material &m, const double &eg, const double &ec);

  /*
    Calculate the integral of the differential cross-section.
    \param m:  material inducing the bremsstrahlung.
    \param eg: energy of the emitted photon.
    \param ec: energy of the charged radiator.
  */
  double isigma(const Material &m, const double &eg, const double &ec);

  // Members.
  double egmin; ///< Minimum allowed photon energy.

private:
  // Members.
  TRandom *rnd; ///< Random number generator.
};

//=============================================================================
/*
  Base class for sub-detectors.
*/
class SubDetector {
public:

  /*
    Base constructor.
    \param mIn: sub-detector material, see Material.
    \param wIn: width of the sub-detector material in mm.
  */
  SubDetector(string mIn, double wIn);

  /*
    Base destructor.
  */
  virtual ~SubDetector();

  /*
    Calculate the z-distance between a particle position and its
    intersection with the sub-detector. If negative, no intersection
    was found. The intersection point is cached.
    \param p: particle to calculate the distance.
   */
  double distance(const Particle &p);
  
  /*
    Intersect method which should be defined for each derived
    sub-detector class.
    \param p: particle to calculate the intersect.
  */
  virtual TLorentzVector intersect(const Particle &p) = 0;

  /*
    Update a particle, e.g. add hits.
  */
  virtual void update(Particle &p);

  // Members.
  TLorentzVector x; ///< Position of last calculated intersect.
  Material       m; ///< Sub-detector material.
  double         w; ///< Width of the sub-detector material in mm.
};

//=============================================================================
/*
  Class to simulate the VELO foils.
*/
class SubDetectorFoils : public SubDetector {
public:
  /*
    Constructor.
    \param r: run configuration to use for the material.
  */
  SubDetectorFoils(string r);

  // Intersect method.
  TLorentzVector intersect(const Particle &p);

private:
  // Members.
  FoilMaterial d;
};

//=============================================================================
/*
  Class to simulate the VELO modules.
*/
class SubDetectorModules : public SubDetector {
public:
  /*
    Constructor.
    \param r: run configuration to use for the material.
  */
  SubDetectorModules(string r);

  // Intersect method.
  TLorentzVector intersect(const Particle &p);

  // Update method.
  void update(Particle &p);

private:
  // Members.
  ModuleMaterial d;
};

//=============================================================================
/*
  Class to simulate the detector.
*/
class Detector {
public:

  /*
    Constructor.
    \param r: run configuration to use for the material.
  */
  Detector(string r);

  /*
    Destructor.
  */
  ~Detector();

  /*
    Transport a particle through the detector.
    \param e: the event record.
    \param i: index of the particle to transport.
  */
  void transport(Event &e, const int &i);

  /*
    Transport all the particles in an event record through the
    detector.
    \param e: the event record to transport.
  */
  void transport(Event &e);

private:
  // Members.
  vector<SubDetector*> sds; ///< Sub-detectors.
  vector<Process*>     pps; ///< Physics processes to apply during transport.
};

#endif // FASTSIM_H
