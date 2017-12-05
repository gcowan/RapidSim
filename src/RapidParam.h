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
            IP,       //IP to own PV
            SIGMAIP,  //Error on IP to own PV
			ETA,      //Pseudorapidity
			PHI,      //Azimuthal angle
			RAPIDITY, //Rapidity
			GAMMA,    //Relitivistic gamma
			BETA,     //velocity
			THETA,    //angle between particles
			COSTHETA, //cosine of angle between particles
			MCORR,    //corrected mass
			UNKNOWN   //unused

		};

		static RapidParam::ParamType typeFromString(TString str);

		RapidParam(TString name, ParamType type, const std::vector<RapidParticle*>& particles, bool truth)
			: name_(name), type_(type), truth_(truth), particles_(particles),
			  minVal_(0.), maxVal_(0.) {setDefaultMinMax();}

		RapidParam(TString name, ParamType type, RapidParticle* part, bool truth)
			: name_(name), type_(type),
			  truth_(truth), minVal_(0.), maxVal_(0.) 
			{particles_.push_back(part); setDefaultMinMax();}

		RapidParam(ParamType type, bool truth)
			: name_(""), type_(type),
			  truth_(truth), minVal_(0.), maxVal_(0.)
			{setDefaultMinMax();}

		~RapidParam() {}

		double eval();//TODO make virtual and give a warning in the base class
		double eval(const TLorentzVector& mom, std::pair<double,double> ip);
		double eval(const TLorentzVector& mom) {return eval(mom,std::pair<double,double>(0.,0.));}

		TString name();
		TString typeName();
		bool truth() { return truth_; }
		double min() { return minVal_; }//TODO make virtual and give a warning in the base class
		double max() { return maxVal_; }//TODO make virtual and give a warning in the base class

		void getMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max);
		void getMinMax(double& min, double& max, bool recalculate=false);

		//convenience functions for 1, 2 or 3-body parameters
		void getMinMax(RapidParticle* part, double& min, double& max);
		void getMinMax(RapidParticle* partA, RapidParticle* partB, double& min, double& max);
		void getMinMax(RapidParticle* partA, RapidParticle* partB, RapidParticle* partC, double& min, double& max);

	private:
		double evalCorrectedMass();
		double evalTheta();

		void setDefaultMinMax() {setDefaultMinMax(particles_,minVal_,maxVal_);}//TODO remove/make virtual and give a warning in the base class

		void setDefaultMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max);
		void setMassMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max);

		TString name_;
		ParamType type_;
		bool truth_;

		//the following are only used for specific parameters
		////TODO refactor into an inherited class RapidParamSpecific
		std::vector<RapidParticle*> particles_;
		double minVal_;
		double maxVal_;

		TLorentzVector mom_;
};

#endif
