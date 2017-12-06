#include "fastsim.h"

//=============================================================================
// Particle class.

// Constructor.
Particle::Particle(const double &px, const double &py, const double &pz,
		   const double &e, int idIn) : 
  TLorentzVector(px, py, pz, e), id(idIn), nvelo(0) {;}
  
// Print.
ostream &operator<<(ostream &s, const Particle &p) {
  s << scientific << setprecision(2) << setw(10) << p.id << setw(10) << p[0] 
    << setw(10) << p[1] << setw(10) << p[2] << setw(10) << p[3] 
    << " (" << setw(9) << p.M() << ")";
  return s;
}

//=============================================================================
// Event class.

// Print.
ostream &operator<<(ostream &s, const Event &e) {
  s << "---------------------------------------------------------------------\n"
    << setw(10) << "id" << setw(10) << "px" << setw(10) << "py" 
    << setw(10) << "pz" << setw(10) << "e" << setw(13) << "(m)\n";
  for (int i = 0; i < (int)e.size(); ++i) s << e[i] << "\n";
  return s;
}

//=============================================================================
// Material class.

// Constructor.
Material::Material(string name) {
  if (name == "Si") {
    z  = 14;
    na = 4.99605e19;
  } else if (name == "Al") {
    z  = 13;
    na = 6.02404e19;
  }
  double a2(pow(alpha*z, 2));
  fi = log(1194.0/pow(z, 2.0/3.0));
  fe = log(184.15/pow(z, 1.0/3.0));
  fc = a2*(0.20206 + 1/(1 + a2) - a2*(0.0396 - 0.0083*a2 + 0.002*a2*a2));
}

//=============================================================================
// ProcessBrem class.

// Constructor.
ProcessBrem::ProcessBrem() : egmin(0.01), rnd(gRandom) {;}

// Add bremsstrahlung.
void ProcessBrem::process(Event &e, const int &i, const Material &m, 
			  const double &w) {
  if (abs(e[i].id) != 11) return;
  double x(0.), eg, tg, pg, ec(e[i][3]), mc(e[i].M());
  TVector3 uc(e[i].Vect().Unit());
  for (int n = 0; n < 1000 && ec - mc > egmin; ++n) {
    x += -log(rnd->Rndm())/(m.na*(isigma(m, ec-mc, ec)-isigma(m, egmin, ec)));
    if (x > w) break;
    eg = energy(m, ec - mc);
    tg = theta(mc, ec);
    pg = phi();
    ec = ec - eg;
    e.push_back(Particle(sin(tg)*cos(pg), sin(tg)*sin(pg), cos(tg), 1));
    e.back() *= eg;
    e.back().id = 22;
    e.back().RotateUz(uc);
    e[i] -= e.back();
    e[i] *= sqrt(ec*ec - mc*mc)/e[i].P();
    e[i][3] = ec;
    uc = e[i].Vect().Unit();
  }
}

// Sample emitted photon energy.
double ProcessBrem::energy(const Material &m, const double &ec) {
  double xs0(dsigma(m, egmin, ec)), xs1(dsigma(m, ec, ec)),
    p0(xs0), p1(ec/log(xs0/xs1)), eg(0);
  for (int n = 0; n < 1000; ++n) {
    eg = rnd->Exp(p1) + egmin;
    if (eg > ec) continue;
    if (rnd->Rndm() < dsigma(m, eg, ec)/(p0*exp(-(eg - egmin)/p1))) break;
  }
  return eg > ec ? ec : eg;
}

// Sample emitted photon theta.
double ProcessBrem::theta(const double &mc, const double &ec) {
  double u(4 - 8*rnd->Rndm()),
    a(mc*sqrt((ec*(ec + 2*mc))/(mc*mc))/(ec + mc)),
    b(pow((abs(u) + sqrt(4 + u*u))/2, 1.0/3.0));
  b = (u < 0 ? -1 : 1)*(1/b - b);
  return acos((a + b)/(1 + a*b));
}

// Sample emitted photon phi.
double ProcessBrem::phi() {return 2*PI*rnd->Rndm();}

// Differential cross-section.
double ProcessBrem::dsigma(const Material &m, const double &eg, 
			   const double &ec) {
  return (4*m.alpha*m.re*m.re*m.z*(ec*(ec - eg)*(1 + m.z) 
    + 3*(4*ec*ec - 4*ec*eg + 3*eg*eg)*(m.fi + (m.fe - m.fc)*m.z)))
    /(9.*ec*ec*eg);
}

// Integral of cross-section.
double ProcessBrem::isigma(const Material &m, const double &eg, 
			   const double &ec) {
  return (2*m.alpha*m.re*m.re*m.z*(eg*(9*eg*(m.fi - m.fc*m.z + m.fe*m.z) 
    - 2*ec*(1 + 12*m.fi + m.z - 12*m.fc*m.z + 12*m.fe*m.z)) + 2*ec*ec
    *(1 + 12*m.fi + m.z - 12*m.fc*m.z + 12*m.fe*m.z)*log(eg)))
    /(9.*ec*ec);
}

//=============================================================================
// SubDetector class.

// Constructor.
SubDetector::SubDetector(string mIn, double wIn) : m(mIn), w(wIn) {;}

// Destructor.
SubDetector::~SubDetector() {;}

// Calculate the intersection z-distance.
double SubDetector::distance(const Particle &p) {
  x = intersect(p);
  return x[3] == 0 ? -1 : abs(x[2] - p.x[2]);
}
  
// Update a particle.
void SubDetector::update(Particle &p) {++p.nvelo;--p.nvelo;}

//=============================================================================
// SubDetectorFoils class.

// Constructor.
SubDetectorFoils::SubDetectorFoils(string r) : SubDetector("Al", 0.3), d(r) {;}


// Intersect method.
TLorentzVector SubDetectorFoils::intersect(const Particle &p) {
  return d.intersect(p.x[0], p.x[1], p.x[2], p[0], p[1], p[2]);
}

//=============================================================================
// SubDetectorModules class.

// Constructor.
SubDetectorModules::SubDetectorModules(string r) : SubDetector("Si", 0.3), d(r) {;}

// Intersect method.
TLorentzVector SubDetectorModules::intersect(const Particle &p) {
  return d.intersect(p.x[0], p.x[1], p.x[2], p[0], p[1], p[2]);
}

// Update method.
void SubDetectorModules::update(Particle &p) {++p.nvelo;}

//=============================================================================
// Detector class.

// Constructor.
Detector::Detector(string r) {
  sds.push_back(new SubDetectorFoils(r));
  sds.push_back(new SubDetectorModules(r));
  pps.push_back(new ProcessBrem());
}

// Destructor.
Detector::~Detector() {
  for (int s = 0; s < (int)sds.size(); ++s) delete sds[s];
  for (int p = 0; p < (int)pps.size(); ++p) delete pps[p];
}

// Transport a particle.
void Detector::transport(Event &e, const int &i) {
  for (int n = 0; n < 100; ++n) {
    double smin(-1), zmin(numeric_limits<double>::infinity()), z;
    for (int s = 0; s < (int)sds.size(); ++s) {
      if (n != 0 && sds[s]->x[3] == 0) continue;
      z = sds[s]->distance(e[i]);
      if (sds[s]->x[3] != 0 && z < zmin) {smin = s; zmin = z;}
    }
    if (smin == -1) break;
    sds[smin]->update(e[i]);
    e[i].x = sds[smin]->x;
    for (int p = 0; p < (int)pps.size(); ++p)
      pps[p]->process(e, i, sds[smin]->m, sds[smin]->w);
    e[i].x.SetVect((1 + 0.1/e[i].x.P())*e[i].x.Vect());
  }
}

// Transport an event.
void Detector::transport(Event &e) {
  for (int i = 0; i < (int)e.size(); ++i) transport(e, i);
}
