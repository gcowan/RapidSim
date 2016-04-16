#ifndef RAPIDPARTICLEDATA_H
#define RAPIDPARTICLEDATA_H

#include <map>
#include <set>

#include "TString.h"

#include "RooRealVar.h"

class RooRelBreitWigner;
class RooGounarisSakurai;
class RapidParticle;

class RapidParticleData {
	public:
		enum ResLineShape {
			RelBW,
			GS
		};

		static RapidParticleData* getInstance();

		void loadData(TString file);

		double getMass(int id);
		double getWidth(int id);
		double getSpin(int id);
		double getCharge(int id);
		TString getName(int id);
		TString getSanitisedName(int id);
		
		int pdgCode(TString name);

		RapidParticle* makeParticle(int id, RapidParticle* mother);
		RapidParticle* makeParticle(TString name, RapidParticle* mother);

		void setupMass(RapidParticle* part);

	private:
		static RapidParticleData* instance_;

		RapidParticleData() {}
		~RapidParticleData() {}

		//copy constructor and copy assignment operator not implemented
		RapidParticleData( const RapidParticleData& other );
		RapidParticleData& operator=( const RapidParticleData& other );

		void addEntry(int id, TString name, double mass, double width, double spin, double charge, TString lineshape);
		TString sanitiseName(TString name);
		TString makeUniqName(TString name);

		RooRelBreitWigner* makeRelBW(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name); 
		RooGounarisSakurai* makeGS(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name); 

		std::map<int, double> idToMass_;
		std::map<int, double> idToWidth_;
		std::map<int, double> idToSpin_;
		std::map<int, double> idToCharge_;
		std::map<int, TString> idToName_;
		std::map<int, ResLineShape> idToShape_;
		std::map<TString, int> nameToId_;

		//register of particle names that have been used
		std::set<TString> usedNames_;
};

#endif
