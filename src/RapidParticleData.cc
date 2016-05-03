#include "RapidParticleData.h"

#include <fstream>
#include <iostream>

#include "RooRelBreitWigner.h"
#include "RooGounarisSakurai.h"

#include "RapidParticle.h"

RapidParticleData* RapidParticleData::instance_=0;

RapidParticleData* RapidParticleData::getInstance() {
	if(!instance_) {
		instance_ = new RapidParticleData();
		instance_->loadData("../config/particles.dat");
	}
	return instance_;
}

void RapidParticleData::loadData(TString file) {
	std::ifstream fin;
	fin.open(file, std::ifstream::in);
	TString buffer;
	buffer.ReadLine(fin);//ignore title line

	int id;
	TString part, anti;
	double mass, width, spin, charge;
	TString lineshape;

	while(fin.good()) {
		buffer.ReadToken(fin);
		id = buffer.Atoi();
		part.ReadToken(fin);
		anti.ReadToken(fin);
		buffer.ReadToken(fin);
		mass = buffer.Atof();
		buffer.ReadToken(fin);
		width = buffer.Atof();
		buffer.ReadToken(fin);
		spin = buffer.Atof();
		buffer.ReadToken(fin);
		charge = buffer.Atof();
		lineshape.ReadToken(fin);

		if(!fin.good()) break;
		addEntry(id,part,mass,width,spin,charge,lineshape);
		if(anti != "---") addEntry(-1*id,anti,mass,width,spin,-1.*charge,lineshape);

	}
	fin.close();
}

double RapidParticleData::getMass(int id) {
	if(idToMass_.count(id)) {
		return idToMass_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getMass : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getWidth(int id) {
	if(idToWidth_.count(id)) {
		return idToWidth_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getWidth : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getSpin(int id) {
	if(idToSpin_.count(id)) {
		return idToSpin_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getSpin : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getCharge(int id) {
	if(idToCharge_.count(id)) {
		return idToCharge_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getCharge : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

TString RapidParticleData::getName(int id) {
	if(idToName_.count(id)) {
		return idToName_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getName : Unknown particle ID " << id << std::endl;
		return "";
	}
}

TString RapidParticleData::getSanitisedName(int id) {
	TString name = getName(id);
	return sanitiseName(name);
}
		
int RapidParticleData::pdgCode(TString name) {
	if(nameToId_.count(name)) {
		return nameToId_[name];
	} else {
		std::cout << "WARNING in RapidParticleData::getId : Unknown particle name " << name << std::endl;
		return 0.;
	}
}

RapidParticle* RapidParticleData::makeParticle(int id, RapidParticle* mother) {
	
	double mass = getMass(id);
	double charge = getCharge(id);
	TString name = getName(id);
	name = makeUniqName(name);

	return new RapidParticle(id, name, mass, charge, mother);
}

RapidParticle* RapidParticleData::makeParticle(TString name, RapidParticle* mother) {

	int id = pdgCode(name);
	double mass = getMass(id);
	double charge = getCharge(id);
	name = makeUniqName(name);

	return new RapidParticle(id, name, mass, charge, mother);
}

void RapidParticleData::setupMass(RapidParticle* part) {
	unsigned int id = TMath::Abs(part->id());

	//check if this is a resonance that we know how to make
	if(!idToShape_.count(id)) return;
	ResLineShape shape = idToShape_[id];
	if(part->nDaughters() != 2) return;

	unsigned int idA = TMath::Abs(part->daughter(0)->id());
	unsigned int idB = TMath::Abs(part->daughter(1)->id());

	double mA = part->daughter(0)->mass();
	double mB = part->daughter(1)->mass();

	if(idA>idB) {
		int tmpId = idB;
		idB = idA;
		idA = tmpId;
		double tmpMass = mB;
		mB = mA;
		mA = tmpMass;
	}

	TString varName = "m";
	varName += part->name();

	TString name = part->name();
	double mass = getMass(id);
	double width = getWidth(id);
	double spin = getSpin(id);

	if(width<narrowWidth_) {
		std::cout << "INFO in RapidParticleData::setupMass : resonance " << name << " is narrow (" << width << " GeV) and will generated using a fixed mass." << std::endl;
		return;
	}

	double mmin = mass - 100.*width;
	double mmax = mass + 100.*width;

	RooRealVar m(varName,varName,mmin,mmax);
	RooAbsPdf* pdf(0);
	switch(shape) {
		case RapidParticleData::GS:
			pdf = makeGS(m, mass, width, spin, mA, mB, name);
			break;
		default:
			std::cout << "WARNING in RapidParticleData::setupMass : unknown lineshape for " << name << "." << std::endl
			          << "                                        : using a relativistic Breit-Wigner." << std::endl;
		case RapidParticleData::RelBW:
			pdf = makeRelBW(m, mass, width, spin, mA, mB, name);
	}
	RooDataSet* massdata = pdf->generate(RooArgSet(m),100000);
	massdata->getRange(m,mmin,mmax);
	part->setMassShape(massdata,mmin,mmax,varName);
}

void RapidParticleData::addEntry(int id, TString name, double mass, double width, double spin, double charge, TString lineshape) {
	idToMass_[id] = mass;
	idToWidth_[id] = width;
	idToSpin_[id] = spin;
	idToCharge_[id] = charge;
	idToName_[id] = name;
	nameToId_[name] = id;
	if(lineshape=="RBW") idToShape_[id] = RapidParticleData::RelBW;
	else if(lineshape=="GS") idToShape_[id] = RapidParticleData::GS;
}

TString RapidParticleData::sanitiseName(TString name) {
	name = name.ReplaceAll("+","p");
	name = name.ReplaceAll("-","m");
	name = name.ReplaceAll("*","st");
	name = name.ReplaceAll("(","_");
	name = name.ReplaceAll(")","_");
	name = name.ReplaceAll("[","_");
	name = name.ReplaceAll("]","_");
	name = name.ReplaceAll("<","_");
	name = name.ReplaceAll(">","_");
	name = name.ReplaceAll("{","_");
	name = name.ReplaceAll("}","_");
	name = name.ReplaceAll(" ","_");
	name = name.ReplaceAll("$","");
	name = name.ReplaceAll("%","");
	name = name.ReplaceAll("&","");
	name = name.ReplaceAll("/","");
	name = name.ReplaceAll(":","");
	name = name.ReplaceAll(";","");
	name = name.ReplaceAll("=","");
	name = name.ReplaceAll("\\","");
	name = name.ReplaceAll("^","");
	name = name.ReplaceAll("|","");
	name = name.ReplaceAll(",","");
	name = name.ReplaceAll(".","");
	name.Remove(TString::kBoth,'_');

	return name;

}

TString RapidParticleData::makeUniqName(TString name) {

	TString uniqName("");

	int i=-1;
	name = sanitiseName(name);

	do {
		++i;
		uniqName = name;
		uniqName+= "_";
		uniqName+= i;
	} while(usedNames_.count(uniqName)>0);

	usedNames_.insert(uniqName);

	return uniqName;

}

RooRelBreitWigner* RapidParticleData::makeRelBW(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name) {

	double barrierFactor=3.;//TODO

	RooRealVar* m0     = new RooRealVar(name+"m0",    name+"m0",     mean);
	RooRealVar* g0     = new RooRealVar(name+"g0",    name+"g0",    gamma);
	RooRealVar* spin   = new RooRealVar(name+"spin",  name+"spin",  thespin);
	RooRealVar* radius = new RooRealVar(name+"radius",name+"radius",barrierFactor); // not used
	RooRealVar* ma     = new RooRealVar(name+"ma",    name+"ma",     m1);
	RooRealVar* mb     = new RooRealVar(name+"mb",    name+"mb",     m2);   
	
	return new RooRelBreitWigner(name,name, m,*m0,*g0,*radius,*ma,*mb,*spin);
}

RooGounarisSakurai* RapidParticleData::makeGS(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name) {

	double barrierFactor=3.;//TODO

	RooRealVar* m0     = new RooRealVar(name+"m0",    name+"m0",     mean);
	RooRealVar* g0     = new RooRealVar(name+"g0",    name+"g0",    gamma);
	RooRealVar* spin   = new RooRealVar(name+"spin",  name+"spin",  thespin);
	RooRealVar* radius = new RooRealVar(name+"radius",name+"radius",barrierFactor); // not used
	RooRealVar* ma     = new RooRealVar(name+"ma",    name+"ma",     m1);
	RooRealVar* mb     = new RooRealVar(name+"mb",    name+"mb",     m2);   
	
	return new RooGounarisSakurai(name,name, m,*m0,*g0,*spin, *radius,*ma,*mb);
}
