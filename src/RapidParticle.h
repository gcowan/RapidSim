#ifndef RAPIDPARTICLE_H
#define RAPIDPARTICLE_H

#include <iostream>
#include <vector>
#include <map>

#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "RooDataSet.h"

#include "functions.h"
#include "RapidMomentumSmear.h"

class RapidParticle {
	public:
		RapidParticle(int id, TString name, RapidParticle* mother)
			: id_(id), name_(name), mass_(0.),
			  mother_(mother), next_(0), stable_(true), invisible_(false)
			{mass_ = getMass(id);}


		~RapidParticle() {}
		
		void addDaughter(RapidParticle* part);
		
		bool generate();
		TLorentzVector& getP() { return p_; }
		TLorentzVector& getPSmeared() { return pSmeared_; }
		
		void smearMomentum();

		int id() { return id_; }
		TString name() { return name_; }
		double mass() { return mass_; }

		unsigned int nDaughters() { return daughters_.size(); }

		RapidParticle* daughter(unsigned int i);
		RapidParticle* next() { return next_; }
		RapidParticle* mother() { return mother_; }

		double* daughterMasses() { return &daughterMasses_[0]; }

		bool stable() { return stable_; }
		bool invisible() { return invisible_; }

		void setId(int id) { id_=id; }
		void setName(TString name) { name_=name; }
		void setMass(double mass) { mass_=mass; }

		void setStable(bool stable=true) { stable_ = stable; }
		void setInvisible(bool invisible=true) { invisible_ = invisible; }
		void setSmearing(RapidMomentumSmear* momSmear) { momSmear_ = momSmear; }

		void setP(TLorentzVector p) { p_ = p; }
		void setPtEtaPhi(double pt, double eta, double phi) { p_.SetPtEtaPhiM(pt,eta,phi,mass_); }

		void print(int index); 

	private:
		int id_;
		TString name_;
		double mass_;
		
		RapidParticle* mother_;
		std::vector<RapidParticle*> daughters_;
		RapidParticle* next_;
		std::vector<double> daughterMasses_;

		TLorentzVector p_;
		TLorentzVector pSmeared_;

		bool stable_;
		bool invisible_;
		
		RapidMomentumSmear* momSmear_;

};
#endif
