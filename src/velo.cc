// Written by Philip Ilten and Mike Williams, 20/06/2017.
// Toolset to check if a secondary vertex is consistent with having
// originated from the VELO material. See README for an overview.
#include "velo.h"

//=============================================================================
// ModuleMaterial::Sensor class.

// Default constructor.
ModuleMaterial::Sensor::Sensor() {;}

// Primary constructor.
ModuleMaterial::Sensor::Sensor(TFile &file, string mod, string sns) : 
  z(sns == "0" ? 4 : 7), half(0) {
  TObjString *exp; TVectorD *par;
  par = (TVectorD*)file.Get((mod + "_a_z_par").c_str());
  if (!par || par->GetNrows() < z + 1) return;
  z   = (*par)[z];
  exp = (TObjString*)file.Get((mod + "_" + sns + "_l_fnc").c_str());
  par = (TVectorD*)file.Get((mod + "_" + sns + "_l_par").c_str());  
  if (!exp || !par) return;
  fncl = TF1(mod.c_str(), exp->GetString().Data(), -100, 100);
  for (int idx = 0; idx < par->GetNrows(); ++idx)
    fncl.SetParameter(idx, (*par)[idx]);
  exp = (TObjString*)file.Get((mod + "_" + sns + "_u_fnc").c_str());
  par = (TVectorD*)file.Get((mod + "_" + sns + "_u_par").c_str());  
  if (!exp || !par) return;
  fncu = TF1(mod.c_str(), exp->GetString().Data(), -100, 100);
  for (int idx = 0; idx < par->GetNrows(); ++idx)
      fncu.SetParameter(idx, (*par)[idx]);
  half = fncl.GetNpar() > fncu.GetNpar() ? 1 : -1;
};

// Operator for sorting by the given z-value.
bool ModuleMaterial::Sensor::operator<(const Sensor& sns) const {
  return z < sns.z;}

// Configure the tool.
void ModuleMaterial::Sensor::config(double xc, double yc, bool rel) {
  if (rel) {yc = fncl.GetParameter(0) + yc; xc = fncl.GetParameter(1) + xc;}
  fncl.SetParameter(0, yc); fncl.SetParameter(1, xc);
  fncu.SetParameter(0, yc); fncu.SetParameter(1, xc);
}

// Configure the tool.
void ModuleMaterial::Sensor::config(const TF1 &fnc) {
  if (half > 0) fncl = TF1(fnc);
  else fncu = TF1(fnc);
}

// Return if an x,y-point falls within a module sensor.
bool ModuleMaterial::Sensor::inside(double x, double y) {
  if (x < fncl.Eval(y)) return false;
  if (x > fncu.Eval(y)) return false;
  return true;
}

// Return the x-position for the edge of the sensor.
double ModuleMaterial::Sensor::x(double y, int edge) {
  if (edge == 0) return (half > 0 ? fncl : fncu).Eval(y);
  else return (edge == -1 ? fncl : fncu).Eval(y);
}

// Return a parameter for the sensor function.
double ModuleMaterial::Sensor::p(unsigned int par, int edge) {
  if (edge == 0) return (half > 0 ? fncl : fncu).GetParameter(par);
  else return (edge == -1 ? fncl : fncu).GetParameter(par);
}

//=============================================================================
// ModuleMaterial::Integral class.

// Return the probability integral for a given half.
double ModuleMaterial::Integral::integral(int half) {
  vector<double> xy0(2, -sigma), xy1(2, sigma);
  double pre(0), vzmin(vz - abs(sigma*dvz)), vzmax(vz + abs(sigma*dvz));
  double dmin(numeric_limits<double>::infinity());
  geo = 0;
  for (int sns = 0; sns < (int)mzs->size(); ++sns) {
    if (mzs->at(sns).half != half) continue;
    double mz(mzs->at(sns).z);
    if (abs(mz - vz) < dmin) {dmin = abs(vz - mz); geo = &mzs->at(sns);}
    if (vzmax < mz) break;
    if (vzmin > mz) continue;
    pre += ((TMath::Erf(sqrt(0.5)*(vz - (mz - width/2))/dvz) -
	     TMath::Erf(sqrt(0.5)*(vz - (mz + width/2))/dvz))/2);
  }
  if (!geo) return 0;
  else return pre*inr->Integral(&xy0[0], &xy1[0]);
}

// Constructor.
ModuleMaterial::Integral::Integral(vector<Sensor> *mzsIn) : 
  mzs(mzsIn), fnc("", this, -100, 100, -100, 100, 0, ""), wrp(fnc), inr(0) {;}

// Destructor.
ModuleMaterial::Integral::~Integral() {if (inr) delete inr;} 

// Return the probability function for a given half.
double ModuleMaterial::Integral::operator()(double *t, double *) {
  double sx(t[0]), sy(t[1]);
  if (geo->inside(sx*dvx + vx, sy*dvy + vy))
    return TMath::Gaus(sx, 0, 1, true)*TMath::Gaus(sy, 0, 1, true);
  else return 0;
}
    
// Configure the tool.
void ModuleMaterial::Integral::config
(double dvminIn, double sigmaIn, double widthIn, double tolIn, 
 ROOT::Math::IntegrationMultiDim::Type algIn) {
  dvmin = abs(dvminIn); sigma = abs(sigmaIn);
  width = abs(widthIn); tol = abs(tolIn);
  if (alg != algIn || !inr) {
    if (inr) delete inr;
    alg = algIn; inr = new ROOT::Math::IntegratorMultiDim(alg);
    inr->SetFunction(wrp);
  }
  inr->SetRelTolerance(tol);
}

// Return the probability integral.
double ModuleMaterial::Integral::integral
(double vxIn, double vyIn, double vzIn, 
 double dvxIn, double dvyIn, double dvzIn, int halfIn) {
  vx = vxIn; vy = vyIn; vz = vzIn;
  dvx = abs(dvxIn); dvy = abs(dvyIn), dvz = abs(dvzIn);
  dvx = dvx > dvmin ? dvx : dvmin;
  dvy = dvy > dvmin ? dvy : dvmin;
  dvz = dvz > dvmin ? dvz : dvmin;
  if (dvx == 0 || dvy == 0 || dvz == 0) return 0;
  double ingr(0);
  if (halfIn == 0 || halfIn == -1) ingr += integral(-1);
  if (halfIn == 0 || halfIn == -1) ingr += integral( 1);
  return ingr > 1 ? 1 : ingr;
}

//=============================================================================
// ModuleMaterial::Distance class.

// Return the distance for a given half.
void ModuleMaterial::Distance::distance(int half, TLorentzVector &s) {
  geo = 0;
  for (int sns = 0; sns < (int)mzs->size(); ++sns) {
    if (mzs->at(sns).half != half) continue;
    if (abs(mzs->at(sns).z - vz) < s.T()) {
      geo = &mzs->at(sns); 
      s.SetXYZT(0, 0, mzs->at(sns).z, abs(mzs->at(sns).z - vz));
    }
  }
  double x(vx), y(vy), z(s.Z()), ymin, ymax;
  if (geo && !geo->inside(vx, vy)) {
    fnc.Update();
    ymin = (geo->half > 0 ? geo->fncl : geo->fncu).GetParameter(0);
    ymax = (geo->half < 0 ? geo->fncl : geo->fncu).GetParameter(2);
    if (y < ymin) {ymax = -ymax; swap(ymin, ymax);}
    y = fnc.GetMinimumX(ymin, ymax, tol, itr, logt);
    x = geo->x(y, 0);
  }
  s.SetXYZT(x, y, z, 0);
  x = (x - vx)/dvx; y = (y - vy)/dvy; z = (z - vz)/dvz; 
  s.SetT(sqrt(x*x + y*y + z*z));
}
 
// Default constructor.
ModuleMaterial::Distance::Distance(vector<Sensor> *mzsIn) : 
  mzs(mzsIn), fnc("", this, -50, 50, 0, "") {fnc.SetNpx(1000);}
    
// Return the transverse distance.
double ModuleMaterial::Distance::operator()(double *t, double *) {
  double x((geo->x(t[0], 0) - vx)/dvx), y((t[0] - vy)/dvy);
  return sqrt(x*x + y*y);
}

// Configure the tool.
void ModuleMaterial::Distance::config
(double dvminIn, double tolIn, int itrIn, bool logtIn) {
  dvmin = abs(dvminIn); tol = abs(tolIn); itr = abs(itrIn); logt = logtIn;
}

// Return the distance.
TLorentzVector ModuleMaterial::Distance::distance
(double vxIn, double vyIn, double vzIn, double dvxIn, double dvyIn,
 double dvzIn, int halfIn) {
  vx = vxIn; vy = vyIn; vz = vzIn; 
  dvx = abs(dvx); dvy = abs(dvy), dvz = abs(dvz);
  dvx = dvxIn > dvmin ? dvxIn : dvmin;
  dvy = dvyIn > dvmin ? dvyIn : dvmin;
  dvz = dvzIn > dvmin ? dvzIn : dvmin;
  TLorentzVector sl(0, 0, 0, numeric_limits<double>::infinity()),
    su(0, 0, 0, numeric_limits<double>::infinity());
  if (dvx == 0 || dvy == 0 || dvz == 0) return sl;
  if (halfIn == 0 || halfIn == -1) distance(-1, sl);
  if (halfIn == 0 || halfIn ==  1) distance( 1, su);
  return sl.T() < su.T() ? sl : su;
}

//=============================================================================
// ModuleMaterial class.

// Constructor.
ModuleMaterial::ModuleMaterial(string pars) : ingr(&mzs), dist(&mzs) {
  configIntegral(); configDistance();
  TFile file(pars.c_str());
  TList *keys = file.GetListOfKeys();
  for (TObjLink *lnk = keys->FirstLink(); lnk; lnk = lnk->Next()) {
    TKey *key((TKey*)lnk->GetObject());
    string name(key->GetName()), mod, sns;
    if (name.find("module") == string::npos) continue;
    if (name.find("_l_fnc") == string::npos) continue;
    mod = name.substr(0, name.length() - 8);
    sns = name.substr(name.length() - 7, 1);
    mzs.push_back(Sensor(file, mod, sns));
  }
  file.Close();
  sort(mzs.begin(), mzs.end());
}

// Return the index for the sensor nearest a given z-position.
int ModuleMaterial::sensor(double z, int half) {
  int idx(0); double dmin(numeric_limits<double>::infinity());
  for (int sns = 0; sns < (int)mzs.size(); ++sns) {
    if (half != 0 && mzs[sns].half != half) continue;
    if (abs(mzs[sns].z - z) < dmin) {idx = sns; dmin = abs(mzs[sns].z - z);}
  }
  return idx;
}

// Return if an x,y-point falls within a module sensor.
bool ModuleMaterial::inside(unsigned int idx, double x, double y) {
  if (idx > mzs.size()) return false;
  return mzs[idx].inside(x, y);
}

// Return the x-position for the edge of a module sensor.
double ModuleMaterial::x(unsigned int idx, double y, int edge) {
  if (idx > mzs.size()) return 0;
  else return mzs[idx].x(y, edge);
}

// Return the z-position for a module sensor.
double ModuleMaterial::z(unsigned int idx) {
  if (idx > mzs.size()) return 0;
  else return mzs[idx].z;
}

// Return a parameter for a module sensor function.
double ModuleMaterial::p(unsigned int idx, int par, int edge) {
  if (idx > mzs.size()) return 0;
  else return mzs[idx].p(par, edge);
}

// Configure the technical sensor parameters.
void ModuleMaterial::configSensor(double xc, double yc, bool rel, int half) {
  for (int sns = 0; sns < (int)mzs.size(); ++sns) {
    if (half != 0 && mzs[sns].half != half) continue;
    mzs[sns].config(xc, yc, rel);
  }
}

// Configure the technical sensor parameters.
void ModuleMaterial::configSensor(FoilMaterial &foil) {
  for (int sns = 0; sns < (int)mzs.size(); ++sns) {
    vector<FoilMaterial::Segment> *fzs = 
      mzs[sns].half < 0 ? &foil.flzs : &foil.fuzs;
    for (int seg = 0; seg < (int)fzs->size(); ++seg)
      if ((*fzs)[seg].valid(mzs[sns].z)) {
	(*fzs)[seg].x(0, mzs[sns].z); mzs[sns].config((*fzs)[seg].fnc); break;}
  }
}

// Configure the technical probability integral parameters.
void ModuleMaterial::configIntegral
(double dvmin, double sigma, double width, double tol,
 ROOT::Math::IntegrationMultiDim::Type alg) {
  ingr.config(dvmin, sigma, width, tol, alg);
}

// Configure the technical distance parameters.
void ModuleMaterial::configDistance
(double dvmin, double tol, int itr, bool logt) {
  dist.config(dvmin, tol, itr, logt);
}

// Return the probability integral.
double ModuleMaterial::integral(double vx, double vy, double vz,
				double dvx, double dvy, double dvz, int half) {
  return ingr.integral(vx, vy, vz, dvx, dvy, dvz, half);
}

// Return the intersection.
TLorentzVector ModuleMaterial::intersect(double vx, double vy, double vz,
					 double fx, double fy, double fz,
					 int half) {
  TLorentzVector s(0, 0, 0, 0);
  if (fz == 0) return s;
  double f(sqrt(fx*fx + fy*fy + fz*fz)); f = f == 0 ? 1 : f;
  fx /= f; fy /= f; fz /= f;
  for (int sns = 0; sns < (int)mzs.size(); ++sns) {
    if (half != 0 && mzs[sns].half != half) continue;
    double mz(mzs[sns].z), t((mz - vz)/fz);
    if (t > 0 && mzs[sns].inside(vx + fx*t, vy + fy*t) && 
	(s.T() == 0 || t < abs(s.T()))) 
      s.SetXYZT(vx + fx*t, vy + fy*t, mz, mzs[sns].half*t);
  }
  if (s.T() != 0) s.SetT(s.T() > 0 ? 1 : -1);
  return s;
}

// Return the distance.
TLorentzVector ModuleMaterial::distance(double vx, double vy, double vz,
					double dvx, double dvy, double dvz,
					int half) {
  return dist.distance(vx, vy, vz, dvx, dvy, dvz, half);
}

//=============================================================================
// FoilMaterial::Segment class.

// Default constructor.
FoilMaterial::Segment::Segment() {;}

// Primary constructor.
FoilMaterial::Segment::Segment(TFile &file, string pre) : zmin(-1), zmax(-1) {
  TObjString *exp((TObjString*)file.Get((pre + "_fnc").c_str()));
  TVectorD *lim((TVectorD*)file.Get((pre + "_lim").c_str())), *par(0);
  TGraph *spl(0);
  if (!exp || !lim || lim->GetNrows() != 2) return;
  zmin = (*lim)[0]; zmax = (*lim)[1];
  fnc = TF1(pre.c_str(), exp->GetString().Data(), -100, 100); pre += "_";
  for (int idx = 0; idx < 6; ++idx) {
    stringstream stridx; stridx << idx; exp = 0; par = 0;
    exp = (TObjString*)file.Get((pre + stridx.str() + "_fnc").c_str());
    par = (TVectorD*)file.Get((pre + stridx.str() + "_par").c_str());
    spl = (TGraph*)file.Get((pre.substr(0, pre.size() - 2) + "a_" 
			     + stridx.str() + "_spl").c_str());
    if (!exp || !par) return;
    fncs.push_back(TF1((pre + stridx.str()).c_str(),
		       exp->GetString().Data(), zmin, zmax));
    if (spl)
      spls.push_back(TSpline3((pre + stridx.str() + "_par").c_str(), spl));
    else 
      spls.push_back(TSpline3());
    mths.push_back(0);
    for (int jdx = 0; jdx < par->GetNrows(); ++jdx)
      fncs[idx].SetParameter(jdx, (*par)[jdx]);
  }
};
    
// Configure the tool.
void FoilMaterial::Segment::config(int par, int mth) {
  if (abs(par) < (int)fncs.size()) mths[par] = mth;
}

// Return if a z-position is in the segment.
bool FoilMaterial::Segment::valid(double &z) {return z >= zmin && z < zmax;}

// Return the x-position for the segment.
double FoilMaterial::Segment::x(double y, double z) {
  for (int idx = 0; idx < (int)fncs.size(); ++idx) {
    if (mths[idx] == 1 && z < spls[idx].GetXmax() && z > spls[idx].GetXmin())
      fnc.SetParameter(idx, spls[idx].Eval(z));
    else fnc.SetParameter(idx, fncs[idx].Eval(z));
  }
  return fnc.Eval(y);
}

// Return a parameter for the segment function.
double FoilMaterial::Segment::p(unsigned int par, double z) {
  if (par >= fncs.size()) return 0;
  if (mths[par] == 1 && z < spls[par].GetXmax() && z > spls[par].GetXmin())
    return spls[par].Eval(z);
  return fncs[par].Eval(z);
}
  
//=============================================================================
// FoilMaterial::Integral class.

// Return the probability integral in one dimension.
double FoilMaterial::Integral::integral
(double v, double dv, double t, double dt) {
  return (TMath::Erf(sqrt(0.5)*(v - (t - dt/2))/dv) -
	  TMath::Erf(sqrt(0.5)*(v - (t + dt/2))/dv))/2;
}
    
// Default constructor.
FoilMaterial::Integral::Integral(FoilMaterial *foilIn) : 
  foil(foilIn), fnc("", this, -100, 100, -100, 100, 0, ""), wrp(fnc), inr(0) {;}

// Destructor.
FoilMaterial::Integral::~Integral() {if (inr) delete inr;}

// Return the probability function for a given half.
double FoilMaterial::Integral::operator()(double *t, double *) {
  double sy(t[0]), sz(t[1]);
  double y(vy + sy*dvy), z(vz + sz*dvz), xf(foil->x(y, z, half));
  return integral(vx, dvx, xf, width)*
    TMath::Gaus(sy, 0, 1, true)*TMath::Gaus(sz, 0, 1, true);
}

// Configure the tool.
void FoilMaterial::Integral::config
(double dvminIn, double sigmaIn, double stepIn, double widthIn, double tolIn,
 ROOT::Math::IntegrationMultiDim::Type algIn) {
  dvmin = abs(dvminIn); sigma = abs(sigmaIn); step = abs(stepIn);
  width = abs(widthIn); tol = abs(tolIn);
  if (alg != algIn || !inr) {
    if (inr) delete inr;
    alg = algIn; inr = new ROOT::Math::IntegratorMultiDim(alg);
    inr->SetFunction(wrp);
  }
  inr->SetRelTolerance(tol);
}

// Return the proibability integral.
double FoilMaterial::Integral::integral(double vxIn, double vyIn, double vzIn,
			      double dvxIn, double dvyIn, double dvzIn,
			      int method) {
  vx = vxIn; vy = vyIn; vz = vzIn;
  dvx = abs(dvxIn); dvy = abs(dvyIn); dvz = abs(dvzIn);
  dvx = dvx > dvmin ? dvx : dvmin; dvy = dvy > dvmin ? dvy : dvmin;
  dvz = dvz > dvmin ? dvz : dvmin;
  double ingr(0);
  if (method == 0) {
    vector<double> yz0(2, -sigma), yz1(2, sigma); 
    for (int z0 = -sigma; z0 < sigma; ++z0) {
      yz0[1] = z0; yz1[1] = z0 + step;
      half = -1; ingr += inr->Integral(&yz0[0], &yz1[0]);
      half =  1; ingr += inr->Integral(&yz0[0], &yz1[0]);
    }
  } else if (method == 1) {
    int yn(2*ceil(sigma*dvy/step)+1), zn(2*ceil(sigma*dvz/step)+1);
    double ymin(vy-(yn-1)/2*step), zmin(vz-(zn-1)/2*step), 
      yds[yn], zd, y, z(zmin);
    for (int yi = 0; yi < yn; ++yi)
      yds[yi] = integral(vy, dvy, ymin + yi*step, step);
    for (int zi = 0; zi < zn; ++zi) {
      y = ymin; zd = integral(vz, dvz, z, step);
      for (int yi = 0; yi < yn && zd != 0; ++yi) {
	if (yds[yi] != 0) {
	  ingr += zd*yds[yi]*(integral(vx, dvx, foil->x(y, z, -1), width) + 
			      integral(vx, dvx, foil->x(y, z,  1), width));
	} y += step;
      } z += step; 
    }
  } else if (method == 2) {
    ingr = TMath::Gaus(foil->x(vy, vz, -1)/dvx, 0, 1, true) + 
      TMath::Gaus(foil->x(vy, vz, 1)/dvx, 0, 1, true);
  }
  return ingr > 1 ? 1 : ingr;
}

//=============================================================================
// FoilMaterial::Intersect class.

// Return the intersection for a given half.
double FoilMaterial::Intersect::intersect(int halfIn, double t) {
  half = halfIn;
  double s = fnc.Eval(t);
  if (s > step/2) return 2e3;
  t = fnc.GetMinimumX(t - step/2, t + step/2, tol, itr, logt);
  s = fnc.Eval(t);
  return (s != s || s > 0.01) ? 2e3 : t;
}

// Check if a point is within valid bounds.
bool FoilMaterial::Intersect::check(double t) {
  double x(vx + fx*t), y(vy + fy*t), z(vz + fz*t);
  if (x < -15 || x > 15) return false;
  if (y < -50 || y > 50) return false;
  if (z < -316 || z > 750) return false; 
  return true;
}

// Default constructor.
FoilMaterial::Intersect::Intersect(FoilMaterial *foilIn) : foil(foilIn) {
  stringstream str; str << this;
  fnc = TF1(str.str().c_str(), this, 0, 1000, 0, "");
  fnc.SetNpx(4);
}
    
// Return the distance from the foil along the flight direction.
double FoilMaterial::Intersect::operator()(double *t, double *) {
  return abs(foil->x(vy+fy*t[0], vz+fz*t[0], half) - (vx+fx*t[0]));
}

// Configure the tool.
void FoilMaterial::Intersect::config
(double stepIn, double tolIn, int itrIn, bool logtIn) {
  step = abs(stepIn); tol = abs(tolIn); itr = abs(itrIn); logt = logtIn;
};

// Return the intersection.
TLorentzVector FoilMaterial::Intersect::intersect
(double vxIn, double vyIn, double vzIn, double fxIn, double fyIn, double fzIn,
 int halfIn) {
  double f(sqrt(fxIn*fxIn + fyIn*fyIn + fzIn*fzIn)), tl(2e3), tu(2e3);
  if (f == 0) return TLorentzVector(0, 0, 0, 0);
  vx = vxIn; vy = vyIn; vz = vzIn; fx = fxIn/f; fy = fyIn/f; fz = fzIn/f;
  for (double t = 0.5; t < 2e3; t += step) {
    if (!check(t)) break;
    if (halfIn == -1 || halfIn == 0) tl = intersect(-1, t);
    if (halfIn ==  1 || halfIn == 0) tu = intersect( 1, t);
    if (tl < tu && tl != 2e3)
      return TLorentzVector(vx + fx*tl, vy + fy*tl, vz + fz*tl, -1);
    else if (tu != 2e3)
      return TLorentzVector(vx + fx*tu, vy + fy*tu, vz + fz*tu, 1);
  }
  return TLorentzVector(0, 0, 0, 0);
}

//=============================================================================
// FoilMaterial::Distance class.

// Return the distance for a given half.
void FoilMaterial::Distance::distance
(int halfIn, int method, int dir, TLorentzVector &s, double &n, double &u) {
  half = halfIn;
  double d, dmin(numeric_limits<double>::infinity()), y(vy), z(vz),
    zmax(vz + (dir != -1)*sigma*(dvz < dvmax ? dvz : dvmax));
  double ymin(foil->p(0, vz, half)), ymax(50);
  double t[2] = {y, z - (dir != 1)*sigma*(dvz < dvmax ? dvz : dvmax) - step};
  if (y < ymin) {ymax = -ymax; ymin += 0.1; swap(ymin, ymax);}
  else ymin -= 0.1;
  min->SetLimitedVariable(0, "y", y, step, ymin, ymax);
  while (t[1] <= zmax) {
    t[1] += step; d = (*this)(t);
    if (d < dmin) {dmin = d; z = t[1];}
    if (method == 2) {
      min->SetVariableValues(t); min->FixVariable(1); min->Minimize();
      u += 1/(min->MinValue() < d ? min->MinValue() : d); ++n;
      min->ReleaseVariable(1);
    }
  }
  min->SetVariableValue(0, vy); min->SetVariableValue(1, vz); min->Minimize();
  if (min->MinValue() < dmin) {
    dmin = min->MinValue(); y = min->X()[0]; z = min->X()[1];}
  min->SetVariableValue(0, y); min->SetVariableValue(1, z); min->Minimize();
  if (min->MinValue() < dmin) {
    dmin = min->MinValue(); y = min->X()[0]; z = min->X()[1];}
  if (dmin < s.T()) s.SetXYZT(foil->x(y, z, half), y, z, dmin);
  if (method == 1) {u += 1/dmin; ++n;}
}
    
// Default constructor.
FoilMaterial::Distance::Distance(FoilMaterial *foilIn) : 
  foil(foilIn), fnc(this, &FoilMaterial::Distance::operator(), 2), min(0) {;}
    
// Return the distance from the foil in uncertainty space.
double FoilMaterial::Distance::operator()(const double *t) {
  double x((vx - foil->x(t[0], t[1], half))/dvx), y((vy - t[0])/dvy),
    z((vz - t[1])/dvz);
  return sqrt(x*x + y*y + z*z);
}

// Configure the tool.
void FoilMaterial::Distance::config
(double dvminIn, double dvmaxIn, double sigmaIn, double stepIn,
 double tolIn, int itrIn, string libIn, string algIn) {
  dvmin = abs(dvminIn); dvmax = abs(dvmaxIn); 
  sigma = abs(sigmaIn); step = abs(stepIn);
  if (libIn != lib || alg != algIn || !min) {
    lib = libIn; alg = algIn;
    min = ROOT::Math::Factory::CreateMinimizer(lib.c_str(), alg.c_str());
  }
  min->SetFunction(fnc); min->SetPrintLevel(0); min->SetValidError(false);
  min->SetMaxIterations(itrIn); min->SetMaxFunctionCalls(itrIn);
  min->SetTolerance(tolIn);
}

// Return the distance.
TLorentzVector FoilMaterial::Distance::distance
(double vxIn, double vyIn, double vzIn, double dvxIn, double dvyIn,
 double dvzIn, int halfIn, int method, int dir) {
  vx = vxIn; vy = vyIn; vz = vzIn;
  dvx = abs(dvxIn); dvy = abs(dvyIn), dvz = abs(dvzIn);
  dvx = dvx > dvmin ? dvx : dvmin;
  dvy = dvy > dvmin ? dvy : dvmin;
  dvz = dvz > dvmin ? dvz : dvmin;
  min->SetLimitedVariable(0, "y", vy, step, vy - sigma*dvy, vy + sigma*dvy);
  min->SetLimitedVariable(1, "z", vz, step,
			  vz - (dir !=  1)*sigma*(dvz < dvmax ? dvz : dvmax),
			  vz + (dir != -1)*sigma*(dvz < dvmax ? dvz : dvmax));
  double n(0), u(0);
  TLorentzVector s(0, 0, 0, numeric_limits<double>::infinity());
  if (dvx == 0 || dvy == 0 || dvz == 0) return s;
  level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;
  if (halfIn == 0 || halfIn == -1) distance(-1, method, dir, s, n, u);
  if (halfIn == 0 || halfIn ==  1) distance( 1, method, dir, s, n, u);
  if (method != 0) s.SetT(n/(2*u));
  gErrorIgnoreLevel = level;
  return s;
}

//=============================================================================
// FoilMaterial class.

// Default constructor.
FoilMaterial::FoilMaterial(string pars) : ingr(this), insc(this), dist(this) {
  configIntegral(); configIntersect(); configDistance();
  vector<string> segs(4, "b"); segs[1] = "c"; segs[2] = "t"; segs[3] = "f";
  TFile file(pars.c_str());
  for (int seg = 0; seg < (int)segs.size(); ++seg) {
    flzs.push_back(Segment(file, "foil_l_" + segs[seg]));
    fuzs.push_back(Segment(file, "foil_u_" + segs[seg]));
  }
  configSegment();
  file.Close();
}
  
// Return the x-position of the foil.
double FoilMaterial::x(double y, double z, int half) {
  vector<Segment> *fzs = half < 0 ? &flzs : &fuzs;
  for (int seg = 0; seg < (int)fzs->size(); ++seg)
    if ((*fzs)[seg].valid(z)) return (*fzs)[seg].x(y, z);
  return 0;
}

// Return a parameter for the foil function.
double FoilMaterial::p(unsigned int par, double z, int half) {
  vector<Segment> *fzs = half < 0 ? &flzs : &fuzs;
  for (int seg = 0; seg < (int)fzs->size(); ++seg)
    if ((*fzs)[seg].valid(z)) return (*fzs)[seg].p(par, z);
  return 0;
}

// Configure the technical segment parameters.
void FoilMaterial::configSegment(int half, string seg, int par, int mth) {
  if (half == 0) {
    configSegment(-1, seg, par, mth); configSegment(1, seg, par, mth); return;
  }
  if (seg == "a") {
    vector<string> segs(4, "b"); segs[1] = "c"; segs[2] = "t"; segs[3] = "f";
    for (int idx = 0; idx < (int)segs.size(); ++idx)
      configSegment(half, segs[idx], par, mth); return;
  }
  if (par < 0) {
    for (int idx = 0; idx < 6; ++idx) configSegment(half, seg, idx, mth);
    return;
  }
  map<string, int> segs; segs["c"] = 1; segs["t"] = 2; segs["f"] = 3;
  (half < 0 ? flzs : fuzs)[segs[seg]].config(par, mth);
}

// Configure the technical probability integral parameters.
void FoilMaterial::configIntegral
(double dvmin, double sigma, double step, double width, double tol,
 ROOT::Math::IntegrationMultiDim::Type alg) {
  ingr.config(dvmin, sigma, step, width, tol, alg);
}


// Configure the technical intersection parameters.
void FoilMaterial::configIntersect
(double step, double tol, int itr, bool logt) {
  insc.config(step, tol, itr, logt);
}

// Configure the technical distance parameters.
void FoilMaterial::configDistance
(double dvmin, double dvmax, double sigma, double step, double tol, int itr, 
 string lib, string alg) {
  dist.config(dvmin, dvmax, sigma, step, tol, itr, lib, alg);
}

// Return the probability integral.
double FoilMaterial::integral
(double vx, double vy, double vz, double dvx, double dvy, double dvz,
 int method) {
  return ingr.integral(vx, vy, vz, dvx, dvy, dvz, method);
}

// Return the intersection.
TLorentzVector FoilMaterial::intersect
(double vx, double vy, double vz, double fx, double fy, double fz, int half) {
  return insc.intersect(vx, vy, vz, fx, fy, fz, half);
}

// Return the distance.
TLorentzVector FoilMaterial::distance
(double vx, double vy, double vz, double dvx, double dvy, double dvz, int half,
 int method, int dir) {
  return dist.distance(vx, vy, vz, dvx, dvy, dvz, half, method, dir);
}

//=============================================================================
// VeloMaterial class.

// Default constructor.
VeloMaterial::VeloMaterial(string pars) : foil(pars), modf(pars), modm(pars) {
  modf.configSensor(foil);
}

// Return the distance.
double VeloMaterial::distance(double vx, double vy, double vz,
			      double dvx, double dvy, double dvz) {
  if (dvz <= 1) dvz = 1;
  if (vz > 300 && dvx < 0.5) dvx = 0.5;
  return 1/(1/modf.distance(vx, vy, vz, dvx, dvy, dvz, -1).T() +
	    1/modf.distance(vx, vy, vz, dvx, dvy, dvz, +1).T() +
	    1/foil.distance(vx, vy, vz, dvx, dvy, dvz, -1, 0, -1).T() +
	    1/foil.distance(vx, vy, vz, dvx, dvy, dvz, -1, 0, +1).T() +
	    1/foil.distance(vx, vy, vz, dvx, dvy, dvz, +1, 0, -1).T() +
	    1/foil.distance(vx, vy, vz, dvx, dvy, dvz, +1, 0, +1).T());
}

// Return true if the flight direction of a vertex intersects a foil tip.
bool VeloMaterial::tip(double vx, double vy, double vz,
		       double fx, double fy, double fz) {
  double mz(modf.z(modf.sensor(vz)));
  return abs(foil.intersect(vx, vy, vz, fx, fy, fz).Z() - mz) < 1 ||
    abs(foil.intersect(vx, vy, vz, -fx, -fy, -fz).Z() - mz) < 1;
}

// Return true if the expected first hit is missed.
bool VeloMaterial::miss(double vx, double vy, double vz,
			double fx, double fy, double fz, double mz) {
  double hz(-numeric_limits<double>::infinity()), hr(0), z(vz);
  while(hr < 8.170){
    TLorentzVector d(modm.intersect(vx + (z - vz)*fx/fz, vy + (z - vz)*fy/fz, 
				    z, fx, fy, fz));
    if(d.T() == 0) break;
    int idx(modm.sensor(d.Z())); 
    double mx(modm.p(idx, 1)), my(modm.p(idx, 0));
    hr = sqrt((d.X() - mx)*(d.X() - mx) + (d.Y() - my)*(d.Y() - my));
    hz = d.Z();
    z  = d.Z() + 1;
  }
  return abs(hz - mz) > 1;
}
