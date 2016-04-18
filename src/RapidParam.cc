#include "RapidParam.h"

#include "TRandom.h"

#include "RapidParticle.h"

double RapidParam::eval() {

	if(type_ == RapidParam::THETA || type_ == RapidParam::COSTHETA) return evalTheta();
	if(type_ == RapidParam::MCORR) return evalCorrectedMass();

	TLorentzVector mom;
	if(truth_) {
		for(unsigned int i=0; i<particles_.size(); ++i) {
			mom += particles_[i]->getP();
		}
	} else {
		for(unsigned int i=0; i<particles_.size(); ++i) {
			mom += particles_[i]->getPSmeared();
		}
	}
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
		case RapidParam::THETA: //dealt with separately above - included to appease compiler
		case RapidParam::COSTHETA: //dealt with separately above - included to appease compiler
		case RapidParam::MCORR: //dealt with separately above - included to appease compiler
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
	double dir = TMath::ACos(cosDir) + gRandom->Gaus(0.,0.01); //TODO smearing parameter should be configurable
	double sinDir = TMath::Sin(dir);
	cosDir = TMath::Cos(dir);

	//get the longitudinal and transverse momentum of the visible daughters wrt the B flight direction
	double pLong  = TMath::Abs(cosDir * momS.P());
	double pTran = TMath::Abs(sinDir * momS.P());

	//invariant masses of the visible daughters and the parent as well as the missing mass
	double mVis2 = momS.M2();
	double mPar2 = momT.M2();
	double mMiss2 = mPar2 - mVis2;

	//the corrected mass
	double mCorr = TMath::Sqrt( mVis2 + pTran*pTran ) + pTran;

	//coefficients of the quadratic equation in pL(invisible)
	double a = 2. * pLong * pLong * mVis2;
	double b = 4*pLong*(2*pTran*pLong - mMiss2);
	double c = 4*pTran*pTran * (pLong*pLong + mPar2) - mMiss2*mMiss2;

	//separate according to whether solutions of pL(invisible) are real or not
	if(b*b - 4.*a*c > 0.) mCorr=-1.*mCorr;

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

TString RapidParam::typeName() {
	switch(type_) {
		case RapidParam::M:
			return "M";
		case RapidParam::M2:
			return "M2";
		case RapidParam::MT:
			return "MT";
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
			return "rapidity";
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
	} else if(str=="rapidity") {
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

void RapidParam::setDefaultMinMax() {
	switch(type_) {
		case RapidParam::M:
		case RapidParam::MCORR:
		case RapidParam::MT:
			minVal_ = 0.;
			maxVal_ = 7.;
			break;
		case RapidParam::M2:
			minVal_ = 0.;
			maxVal_ = 30.;
			break;
		case RapidParam::E:
		case RapidParam::ET:
		case RapidParam::P:
		case RapidParam::PX:
		case RapidParam::PY:
		case RapidParam::PZ:
		case RapidParam::PT:
			minVal_ = 0.;
			maxVal_ = 100.;
			break;
		case RapidParam::ETA:
			minVal_ = -7.;
			maxVal_ =  7.;
			break;
		case RapidParam::PHI:
			minVal_ = -3.5;
			maxVal_ =  3.5;
			break;
		case RapidParam::RAPIDITY:
			minVal_ = 0.;
			maxVal_ = 10.;
			break;
		case RapidParam::GAMMA:
			minVal_ = 1.;
			maxVal_ = 50.;
			break;
		case RapidParam::BETA:
			minVal_ = 0.;
			maxVal_ = 1.;
			break;
		case RapidParam::THETA:
			minVal_ = 0.0;
			maxVal_ = 3.5;
			break;
		case RapidParam::COSTHETA:
			minVal_ = -1.0;
			maxVal_ =  1.0;
			break;
	}
}
