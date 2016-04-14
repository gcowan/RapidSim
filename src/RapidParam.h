#ifndef RAPIDPARAM_H
#define RAPIDPARAM_H

#include <vector>

#include "TLorentzVector.h"

class RapidParticle;

class RapidParam {
	public:
		enum ParamType {
			M,        //Mass
			M2,       //Mass squared
			MT,       //Transverse mass
			E,        //Energy
			ET,       //Transverse energy
			P,        //3-momentum
			PX,       //X-momentum
			PY,       //Y-momentum
			PZ,       //Z-momentum
			PT,       //Transverse momentum
			ETA,      //Pseudorapidity
			PHI,      //Azimuthal angle
			RAPIDITY, //Rapidity
			GAMMA,    //Relitivistic gamma
			BETA,     //velocity
			THETA,    //angle from other particle
			MCORR     //corrected mass

		};

		static RapidParam::ParamType typeFromString(TString str);

		RapidParam(TString name, ParamType type, std::vector<RapidParticle*> particles, bool truth)
			: name_(name), type_(type), particles_(particles),
			  truth_(truth), minVal_(0.), maxVal_(0.) {setDefaultMinMax();}

		RapidParam(TString name, ParamType type, RapidParticle* part, bool truth)
			: name_(name), type_(type),
			  truth_(truth), minVal_(0.), maxVal_(0.) 
			{particles_.push_back(part); setDefaultMinMax();}

		~RapidParam() {}

		double eval();

		TString name() { return name_; }
		TString typeName();
		double min() { return minVal_; }
		double max() { return maxVal_; }

	private:
		double evalCorrectedMass();
		double evalTheta();

		void setDefaultMinMax();

		TString name_;
		ParamType type_;
		std::vector<RapidParticle*> particles_;
		bool truth_;
		double minVal_;
		double maxVal_;
};

#endif
