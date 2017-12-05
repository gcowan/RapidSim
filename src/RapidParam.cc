#include "RapidParam.h"

#include "TRandom.h"

#include "RapidParticle.h"
#include "RapidParticleData.h"

double RapidParam::eval() {

	if(type_ == RapidParam::THETA || type_ == RapidParam::COSTHETA) return evalTheta();
	if(type_ == RapidParam::MCORR) return evalCorrectedMass();

	mom_.SetPxPyPzE(0.,0.,0.,0.);
	if(truth_) {
		for(unsigned int i=0; i<particles_.size(); ++i) {
			mom_ += particles_[i]->getP();
		}
	} else {
		for(unsigned int i=0; i<particles_.size(); ++i) {
			mom_ += particles_[i]->getPSmeared();
		}
	}

	return eval(mom_);
}

double RapidParam::eval(const TLorentzVector& mom, std::pair<double,double> ip) {
	switch(type_) {
		case RapidParam::M:
			return mom.M();
		case RapidParam::M2:
			return mom.M2();
		case RapidParam::MT:
			return mom.Mt();
		case RapidParam::E:
			return mom.E();
		case RapidParam::ET:
			return mom.Et();
		case RapidParam::P:
			return mom.P();
		case RapidParam::PX:
			return mom.Px();
		case RapidParam::PY:
			return mom.Py();
		case RapidParam::PZ:
			return mom.Pz();
		case RapidParam::PT:
			return mom.Pt();
		case RapidParam::ETA:
			return mom.Eta();
		case RapidParam::PHI:
			return mom.Phi();
		case RapidParam::RAPIDITY:
			return mom.Rapidity();
		case RapidParam::GAMMA:
			return mom.Gamma();
		case RapidParam::BETA:
			return mom.Beta();
        case RapidParam::IP:
            return ip.first;
        case RapidParam::SIGMAIP:
            return ip.second;
		case RapidParam::THETA:
		case RapidParam::COSTHETA:
		case RapidParam::MCORR:
			std::cout << "WARNING in RapidParam::eval : parameter type " << type_ << " cannot be used with a single momentum." << std::endl
				  << "                              returning 0." << std::endl;
			break;
		default:
			std::cout << "WARNING in RapidParam::eval : unknown parameter type " << type_ << std::endl
				  << "                              returning 0." << std::endl;

	}

	return 0.;

}

double RapidParam::evalCorrectedMass() {

	TLorentzVector momS, momT;

	//load the true (inc. invisible) and smeared momenta
	for(unsigned int i=0; i<particles_.size(); ++i) {
		momT += particles_[i]->getP();
		momS += particles_[i]->getPSmeared();
	}

	//get the sine and cosine of the angle to the B flight direction with a Gaussian smearing applied
	double cosDir = momT.Vect().Dot(momS.Vect())/(momT.P()*momS.P());
	double dir = TMath::ACos(cosDir) + gRandom->Gaus(0.,0.001); //TODO smearing parameter should be configurable
	double sinDir = TMath::Sin(dir);
	cosDir = TMath::Cos(dir);

	//get the longitudinal and transverse momentum of the visible daughters wrt the B flight direction
	//double pLong  = TMath::Abs(cosDir * momS.P());
	double pTran = TMath::Abs(sinDir * momS.P());

	//invariant masses of the visible daughters and the parent as well as the missing mass
	double mVis2 = momS.M2();
	//double mPar2 = momT.M2();
	//double mMiss2 = mPar2 - mVis2;

	//the corrected mass
	double mCorr = TMath::Sqrt( mVis2 + pTran*pTran ) + pTran;

	return mCorr;
}

double RapidParam::evalTheta() {
	if(particles_.size()<2) {
		return -99.;
	}

	TLorentzVector pA, pB, pBoost;

	if(truth_) {
		pA = particles_[0]->getP();
		pB = particles_[1]->getP();
		for(unsigned int i=2; i<particles_.size(); ++i) {
			pBoost += particles_[i]->getP();
		}
	} else {
		pA = particles_[0]->getPSmeared();
		pB = particles_[1]->getPSmeared();
		for(unsigned int i=2; i<particles_.size(); ++i) {
			pBoost += particles_[i]->getPSmeared();
		}
	}

	pA.Boost(-pBoost.BoostVector());
	pB.Boost(-pBoost.BoostVector());

	if(type_ == RapidParam::THETA) {
		return pA.Angle(pB.Vect());
	} else {
		return pA.Vect().Dot(pB.Vect())/(pA.P()*pB.P());
	}
}

TString RapidParam::name() {
	if(name_.Length()==0) {
		//construct default name
		std::vector<RapidParticle*>::iterator it = particles_.begin();
		for( ; it!= particles_.end(); ++it) {
			name_ += (*it)->name();
			name_ += "_";
		}

		name_ += typeName();

		if(truth_) {
			name_ += "_TRUE";
		}
	}

	return name_;
}

TString RapidParam::typeName() {
	switch(type_) {
		case RapidParam::M:
			return "M";
		case RapidParam::M2:
			return "M2";
		case RapidParam::MT:
			return "MT";
		case RapidParam::IP:
            return "IP";
        case RapidParam::SIGMAIP:
            return "SIGMAIP";
        case RapidParam::E:
			return "E";
		case RapidParam::ET:
			return "ET";
		case RapidParam::P:
			return "P";
		case RapidParam::PX:
			return "PX";
		case RapidParam::PY:
			return "PY";
		case RapidParam::PZ:
			return "PZ";
		case RapidParam::PT:
			return "PT";
		case RapidParam::ETA:
			return "eta";
		case RapidParam::PHI:
			return "phi";
		case RapidParam::RAPIDITY:
			return "y";
		case RapidParam::GAMMA:
			return "gamma";
		case RapidParam::BETA:
			return "beta";
		case RapidParam::THETA:
			return "theta";
		case RapidParam::COSTHETA:
			return "costheta";
		case RapidParam::MCORR:
			return "Mcorr";
		default:
			std::cout << "WARNING in RapidParam::typeName : unknown type " << type_ << "." << std::endl
				  << "                                  returning empty string." << std::endl;
			return "";
	}
}

RapidParam::ParamType RapidParam::typeFromString(TString str) {
	if(str=="M") {
		return RapidParam::M;
	} else if(str=="M2") {
		return RapidParam::M2;
	} else if(str=="MT") {
		return RapidParam::MT;
	} else if(str=="IP") {
        return RapidParam::IP;
    } else if(str=="SIGMAIP") {
        return RapidParam::SIGMAIP;
    } else if(str=="E") {
		return RapidParam::E;
	} else if(str=="ET") {
		return RapidParam::ET;
	} else if(str=="P") {
		return RapidParam::P;
	} else if(str=="PX") {
		return RapidParam::PX;
	} else if(str=="PY") {
		return RapidParam::PY;
	} else if(str=="PZ") {
		return RapidParam::PZ;
	} else if(str=="PT") {
		return RapidParam::PT;
	} else if(str=="eta") {
		return RapidParam::ETA;
	} else if(str=="phi") {
		return RapidParam::PHI;
	} else if(str=="y") {
		return RapidParam::RAPIDITY;
	} else if(str=="gamma") {
		return RapidParam::GAMMA;
	} else if(str=="beta") {
		return RapidParam::BETA;
	} else if(str=="theta") {
		return RapidParam::THETA;
	} else if(str=="costheta") {
		return RapidParam::COSTHETA;
	} else if(str=="Mcorr") {
		return RapidParam::MCORR;
	} else {
		std::cout << "WARNING in RapidParam::typeFromString : unknown type name " << str << "." << std::endl
			  << "                                        returning mass parameter type." << std::endl;
		return RapidParam::M;
	}
}

void RapidParam::getMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max) {
	setDefaultMinMax(parts,min,max);
}

void RapidParam::getMinMax(double& min, double& max, bool recalculate) {
	if(!recalculate) {
		min = minVal_;
		max = maxVal_;
	} else {
		setDefaultMinMax(particles_,min,max);
	}
}

void RapidParam::getMinMax(RapidParticle* part, double& min, double& max) {
	std::vector<RapidParticle*> parts;
	parts.push_back(part);
	setDefaultMinMax(parts,min,max);
}

void RapidParam::getMinMax(RapidParticle* partA, RapidParticle* partB, double& min, double& max) {
	std::vector<RapidParticle*> parts;
	parts.push_back(partA);
	parts.push_back(partB);
	setDefaultMinMax(parts,min,max);
}

void RapidParam::getMinMax(RapidParticle* partA, RapidParticle* partB, RapidParticle* partC, double& min, double& max) {
	std::vector<RapidParticle*> parts;
	parts.push_back(partA);
	parts.push_back(partB);
	parts.push_back(partC);
	setDefaultMinMax(parts,min,max);
}

void RapidParam::setDefaultMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max) {
	switch(type_) {
		case RapidParam::M:
			setMassMinMax(parts,min,max);
			break;
		case RapidParam::MT:
			setMassMinMax(parts,min,max);
			min = 0.;
			break;
		case RapidParam::MCORR:
			setMassMinMax(parts,min,max);
			min = -max;
			break;
		case RapidParam::M2:
			setMassMinMax(parts,min,max);
			min = min*min;
			max = max*max;
			break;
        case RapidParam::IP:
        case RapidParam::SIGMAIP:
		case RapidParam::E:
		case RapidParam::ET:
		case RapidParam::P:
		case RapidParam::PX:
		case RapidParam::PY:
		case RapidParam::PZ:
		case RapidParam::PT:
			min = 0.;
			max = 100.;
			break;
		case RapidParam::ETA:
			min = -7.;
			max =  7.;
			break;
		case RapidParam::PHI:
			min = -3.5;
			max =  3.5;
			break;
		case RapidParam::RAPIDITY:
			min = 0.;
			max = 10.;
			break;
		case RapidParam::GAMMA:
			min = 1.;
			max = 50.;
			break;
		case RapidParam::BETA:
			min = 0.;
			max = 1.;
			break;
		case RapidParam::THETA:
			min = 0.0;
			max = 3.5;
			break;
		case RapidParam::COSTHETA:
			min = -1.0;
			max =  1.0;
			break;
		case RapidParam::UNKNOWN:
		default:
			return;
	}
}

void RapidParam::setMassMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max) {
	if(parts.size()==0) {
		//no particles so mass is 0
		return;
	} else if(parts.size()==1) {
		//single particle is easy
		min = parts[0]->minMass() - 0.1;
		max = parts[0]->maxMass() + 0.1;

		//now check to see if we've used any alternative mass hypotheses and extend range if necessary
		double deltaDown(0.), deltaUp(0.);

		RapidParticleData::getInstance()->getMaxAltHypothesisMassShifts(parts,deltaDown,deltaUp);

		//TODO factor of 2. is arbitrary - there is probably a better way to calculate these
		min+=2.*deltaDown;
		max+=2.*deltaUp;

		if(min<0.) min=0.;
	} else {
		RapidParticleData* rpd = RapidParticleData::getInstance();

		//first check whether any of the particles we've been given have a hierarchical relationship
		if(rpd->checkHierarchy(parts)) {
			std::cout << "WARNING in RapidParam::setMassMinMax : some of the particles used in parameter " << name() << " are repeated or have a hierarchical relationship." << std::endl
			          << "                                       setting default mass range." << std::endl;
			min=0.;
			max=10.;
			return;
		}

		//find the most recent common ancestor of all our particles
		RapidParticle* commonAncestor = rpd->findCommonAncestor(parts);

		//this shouldn't happen
		if(!commonAncestor) {
			std::cout << "WARNING in RapidParam::setMassMinMax : some of the particles used in parameter " << name() << " have no common ancestors." << std::endl
			          << "                                       setting default mass range." << std::endl;
			min=0.;
			max=10.;
			return;
		}

		//get the list of stable particles descended from our particles
		std::vector<RapidParticle*> selStableParts;
		rpd->findStableDaughters(parts,selStableParts);

		//now compare against the list of stable particles from the ancestor to get the list of particles we did not include
		std::vector<RapidParticle*> otherParts;
		rpd->findOtherDaughters(commonAncestor, selStableParts, otherParts);

		//now combine any particles where all of the daughters are in this list as the breakup mass is also lost
		std::vector<RapidParticle*> otherPartsCombined;
		rpd->combineCompleteAncestors(otherParts, otherPartsCombined);

		//minimum mass is the sum of the masses of the included particles
		min = -0.1;
		std::vector<RapidParticle*>::const_iterator cit = parts.begin();
		for( ; cit!=parts.end(); ++cit) {
			min += (*cit)->minMass();
		}

		//maximum mass is the mass of the parent minus the sum of the masses of the missing daughters
		max = commonAncestor->maxMass() + 0.1;
		std::vector<RapidParticle*>::iterator it = otherPartsCombined.begin();
		for( ; it!=otherPartsCombined.end(); ++it) {
			max -= (*it)->minMass();
		}

		//now check to see if we've used any alternative mass hypotheses and extend range if necessary
		double deltaDown(0.), deltaUp(0.);

		RapidParticleData::getInstance()->getMaxAltHypothesisMassShifts(parts,deltaDown,deltaUp);

		//TODO factor of 2. is arbitrary - there is probably a better way to calculate these
		min+=2.*deltaDown;
		max+=2.*deltaUp;

		if(min<0.) min=0.;
	}
}
