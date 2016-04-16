#ifndef RAPIDPARTICLE_H
#define RAPIDPARTICLE_H

#include <vector>

#include "TLorentzVector.h"
#include "TString.h"

#include "RooDataSet.h"

class RapidMomentumSmear;
class RapidParticleData;

class RapidParticle {
	public:
		RapidParticle(int id, TString name, double mass, double charge, RapidParticle* mother)
			: index_(0), id_(id), name_(name), mass_(mass), charge_(charge),
			  mother_(mother), next_(0), stable_(true), invisible_(false),
			  massData_(0), minMass_(mass), maxMass_(mass)
			{}


		~RapidParticle() {}
		
		void addDaughter(RapidParticle* part);
		
		bool generate();
		TLorentzVector& getP() { return p_; }
		TLorentzVector& getPSmeared() { return pSmeared_; }
		
		void smearMomentum();

		int id() { return id_; }
		TString name() { return name_; }
		double mass() { return mass_; }
		double charge() { return charge_; }
		double minMass() { return minMass_; }
		double maxMass() { return maxMass_; }

		unsigned int nDaughters() { return daughters_.size(); }

		RapidParticle* daughter(unsigned int i);
		RapidParticle* next() { return next_; }
		RapidParticle* mother() { return mother_; }

		double* daughterMasses() { return &daughterMasses_[0]; }

		bool stable() { return stable_; }
		bool invisible() { return invisible_; }

		void setId(int id) { id_=id; }
		void setName(TString name) { name_=name; }
		void setMass(double mass);

		void setStable(bool stable=true) { stable_ = stable; }
		void setInvisible(bool invisible=true) { invisible_ = invisible; }
		void setSmearing(RapidMomentumSmear* momSmear) { momSmear_ = momSmear; }

		void setP(TLorentzVector p) { p_ = p; }
		void setPtEtaPhi(double pt, double eta, double phi) { p_.SetPtEtaPhiM(pt,eta,phi,mass_); }

		void print(int index); 

		void setMassShape(RooDataSet* ds, double minMass, double maxMass, TString varName);
		void floatMass();

	private:
		void updateDaughterMass(unsigned int index);

		unsigned int index_;

		int id_;
		TString name_;
		double mass_;
		double charge_;
		
		RapidParticle* mother_;
		std::vector<RapidParticle*> daughters_;
		RapidParticle* next_;
		std::vector<double> daughterMasses_;

		TLorentzVector p_;
		TLorentzVector pSmeared_;

		bool stable_;
		bool invisible_;
		
		RapidMomentumSmear* momSmear_;

		RooDataSet* massData_;
		double minMass_;
		double maxMass_;
		TString varName_;
};
#endif
