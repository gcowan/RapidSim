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
			  mother_(mother), next_(0), invisible_(false), momSmear_(0),
			  massData_(0), minMass_(mass), maxMass_(mass),
			  evtGenModel_("PHSP"),
			  currentHypothesis_(0)
			{setPtEtaPhi(0,0,0);}


		~RapidParticle() {}

		void addDaughter(RapidParticle* part);

		void addMassHypothesis(TString name, double mass);
		void setMassHypothesis(unsigned int i);

		bool generate();
		TLorentzVector& getP() { return p_; }
		TLorentzVector& getPSmeared() { return pSmeared_; }

		void smearMomentum();

		int id() { return id_; }
		TString name() { return name_; }
		double mass() { return mass_; }
		double deltaMass();
		double charge() { return charge_; }
		double minMass() { return minMass_; }
		double maxMass() { return maxMass_; }

		unsigned int nDaughters() { return daughters_.size(); }

		unsigned int nMassHypotheses() { return altMasses_.size()+1; }
		TString massHypothesisName();

		bool hasCharm();
		bool hasBeauty();

		RapidParticle* daughter(unsigned int i);
		RapidParticle* next() { return next_; }
		RapidParticle* mother() { return mother_; }

		double const * daughterMasses() { return &daughterMasses_[0]; }

		bool stable() { return daughters_.empty(); }
		bool invisible() { return invisible_; }

		void setName(TString name) { name_=name; }

		void setInvisible(bool invisible=true) { invisible_ = invisible; }
		void setSmearing(RapidMomentumSmear* momSmear) { momSmear_ = momSmear; }

		void setP(TLorentzVector p) { p_ = p; }
		void setPtEtaPhi(double pt, double eta, double phi) { p_.SetPtEtaPhiM(pt,eta,phi,mass_); }

		void print(int index); 

		void setMassShape(RooDataSet* ds, double minMass, double maxMass, TString varName);
		void floatMass();

		TString evtGenDecayModel() { return evtGenModel_; }
		void setEvtGenDecayModel(TString value) { evtGenModel_ = value; }

	private:
		bool hasFlavour(int flavour);

		void setMass(double mass);
		void updateDaughterMass(unsigned int index);

		void updateMomenta();

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

		bool invisible_;

		RapidMomentumSmear* momSmear_;

		RooDataSet* massData_;
		double minMass_;
		double maxMass_;
		TString varName_;

		TString evtGenModel_;

		//store alternative mass hypotheses
		std::vector<TString> massHypothesisNames_;
		std::vector<double> altMasses_;

		unsigned int currentHypothesis_;
};
#endif
