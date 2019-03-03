#include "RapidDecay.h"

#include <iostream>
#include <vector>

#include "TMath.h"
#include "TRandom.h"

#include "RapidExternalEvtGen.h"
#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParam.h"
#include "RapidParticle.h"
#include "RapidParticleData.h"
#include "RapidBeamData.h"
#include "RapidVertex.h"

void RapidDecay::setParentKinematics(TH1* ptHisto, TH1* etaHisto) {
	std::cout << "INFO in RapidDecay::setParentKinematics : setting kinematics of the parent." << std::endl;
	ptHisto_=ptHisto;
	etaHisto_=etaHisto;
}

void RapidDecay::setAcceptRejectHist(TH1* histo, RapidParam* param) {
	accRejParameterX_ = param;
	accRejParameterY_ = 0;
	accRejHisto_      = histo;

	//correct the histogram to account for the the phasespace distribution
	TH1* denom = generateAccRejDenominator1D();
	accRejHisto_->Divide(denom);
	delete denom;
}

void RapidDecay::setAcceptRejectHist(TH1* histo, RapidParam* paramX, RapidParam* paramY) {
	accRejParameterX_ = paramX;
	accRejParameterY_ = paramY;
	accRejHisto_      = histo;

	//correct the histogram to account for the the phasespace distribution
	TH2* denom = generateAccRejDenominator2D();
	accRejHisto_->Divide(denom);
	delete denom;
}

void RapidDecay::setExternal(RapidExternalGenerator* external) {
	external_ = external;
}

bool RapidDecay::checkDecay() {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		if(part->nDaughters()>0) {
			if(!decay_.SetDecay(part->getP(), part->nDaughters(), part->daughterMasses())) {
				std::cout << "ERROR in RapidDecay::checkDecay : decay of " << part->name() << " is kinematically forbidden." << std::endl;
				return false;
			}
		}
	}
	return true;
}

bool RapidDecay::generate(bool genpar) {
	//keep resonance masses and parent kinematics independent of the accept/reject decision
	//these will only be biased if the function is very inefficient for certain values
	//however, one should not use an a/r function the is highly correlated to these variables
	floatMasses();
	if (genpar) genParent();

	bool decayed(false);
	if(external_) {
		decayed = external_->decay(parts_);
	}

	if(!decayed) {
		if(accRejHisto_) {
			if(!genDecayAccRej()) return false;

		} else {
			if(!genDecay()) return false;
		}
	}

	smearMomenta();

	return true;
}

void RapidDecay::smearMomenta() {
	//run backwards so that we reach the daughters first
	for(int i=parts_.size()-1; i>=0; --i) {//don't change to unsigned - needs to hit -1 to break loop
		parts_[i]->smearMomentum();
	}

}

void RapidDecay::setup() {
	setupMasses();

	std::cout << "INFO in RapidDecay::setup : Particle summary follows:" << std::endl;
	printf("index\tlabel\t\t   ID\t\tmass (GeV/c^2)\tparent\t\t# children\tchildren\n");
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->print(i);
	}
}

void RapidDecay::floatMasses() {
	for(unsigned int i=0; i<parts_.size(); ++i) {
		parts_[i]->floatMass();
	}
}

void RapidDecay::setupMasses() {
	RapidParticleData* particleData = RapidParticleData::getInstance();

	for(unsigned int i=0; i<parts_.size(); ++i) {
		particleData->setupMass(parts_[i]);
	}
}

bool RapidDecay::runAcceptReject() {
	if(accRejParameterY_) return runAcceptReject2D();
	else return runAcceptReject1D();
}

bool RapidDecay::runAcceptReject1D() {
	double val = accRejParameterX_->eval();
	int bin = accRejHisto_->FindBin(val);

	double score(0.);
	if(!accRejHisto_->IsBinOverflow(bin) && !accRejHisto_->IsBinUnderflow(bin)) {
		score = accRejHisto_->Interpolate(val);
	}
	double max = accRejHisto_->GetMaximum();
	if(score > gRandom->Uniform(max)) return true;
	return false;
}

bool RapidDecay::runAcceptReject2D() {
	double valX = accRejParameterX_->eval();
	double valY = accRejParameterY_->eval();
	int bin = accRejHisto_->FindBin(valX,valY);

	double score(0.);
	if(!accRejHisto_->IsBinOverflow(bin) && !accRejHisto_->IsBinUnderflow(bin)) {
		score = accRejHisto_->Interpolate(valX,valY);
	}
	double max = accRejHisto_->GetMaximum();
	if(score > gRandom->Uniform(max)) return true;
	return false;
}

TH1* RapidDecay::generateAccRejDenominator1D() {
	TH1* denomHisto = dynamic_cast<TH1*>(accRejHisto_->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in RapidDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay(true)) continue;
		denomHisto->Fill(accRejParameterX_->eval());
	}
	return denomHisto;
}

TH2* RapidDecay::generateAccRejDenominator2D() {
	TH2* denomHisto = dynamic_cast<TH2*>(accRejHisto_->Clone("denom"));
	denomHisto->Reset();

	std::cout << "INFO in RapidDecay::generateAccRejDenominator : generating 1M decays to remove the \"phasespace\" distribution..." << std::endl;
	for(int i=0; i<1000000; ++i) {
		floatMasses();
		genParent();
		if(!genDecay(true)) continue;
		denomHisto->Fill(accRejParameterX_->eval(), accRejParameterY_->eval());
	}
	return denomHisto;
}

void RapidDecay::genParent() {
	double pt(0), eta(0), phi(gRandom->Uniform(0,2*TMath::Pi()));
	if(ptHisto_)   pt = ptHisto_->GetRandom();
	if(etaHisto_) eta = etaHisto_->GetRandom();
	parts_[0]->setPtEtaPhi(pt,eta,phi);

	//Now the pileup vertices
	RapidBeamData* beam = RapidBeamData::getInstance();
	unsigned int numpileup_ = gRandom->Poisson(beam->getPileup());
	double sigmapvxy_ = beam->getSigmaXY();
	double sigmapvz_  = beam->getSigmaZ();

	pileuppvs_.clear();
	for(unsigned int i=0; i<numpileup_; ++i) {
		RapidVertex vtx(gRandom->Gaus(0,sigmapvxy_),gRandom->Gaus(0,sigmapvxy_),gRandom->Gaus(0,sigmapvz_));
		pileuppvs_.push_back(vtx);
	}
}

bool RapidDecay::genDecay(bool acceptAny) {
	//The origin vertex of the signal is always 0,0,0
	RapidVertex * signalpv = parts_[0]->getOriginVertex();
	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];
		if(part->nDaughters()>0) {
			// check decay kinematics valid
			if(!decay_.SetDecay(part->getP(), part->nDaughters(), part->daughterMasses())) {
				if(!suppressKinematicWarning_) {
					std::cout << "WARNING in RapidDecay::genDecay : decay of " << part->name() << " is kinematically forbidden for some events due to resonance mass shapes." << std::endl
						  << "                                  these events will not be generated." << std::endl
						  << "                                  further warnings will be suppressed." << std::endl;
					suppressKinematicWarning_ = true;
				}
				return false;
			}

			// make an event
			if(acceptAny) {
				decay_.Generate();
			} else {
				int nGen(0);
				bool accept(false);
				while (nGen < maxgen_ && accept == false){
					accept = decay_.Generate() > gRandom->Uniform();
					++nGen;
				} // while

				if(!accept) {
					if(!suppressAttemptsWarning_) {
						std::cout << "WARNING in RapidDecay::genDecay : rejected all " << maxgen_ << " attempts to decay " << part->name() << "." << std::endl
							  << "                                  this event will not be generated." << std::endl
							  << "                                  further warnings will be suppressed." << std::endl;
						suppressAttemptsWarning_ = true;
					}
					return false;
				}
			}
			// Now generate the decay vertex for long-lived particles
			// First set the origin vertex to be the PV for the head of the chain
			// in all other cases, the origin vertex will already be set in the loop below
			if (part->ctau()>0) {
				double dist = part->getP().P()*gRandom->Exp(part->ctau())/part->mass();
				double dvx  = part->getOriginVertex()->getVertex(true).X() + part->getP().Vect().Unit().X()*dist;
				double dvy  = part->getOriginVertex()->getVertex(true).Y() + part->getP().Vect().Unit().Y()*dist;
				double dvz  = part->getOriginVertex()->getVertex(true).Z() + part->getP().Vect().Unit().Z()*dist;
				part->getDecayVertex()->setXYZ(dvx,dvy,dvz,part->getDecayVertexSmear());
			}

			int j=0;
			for(RapidParticle* jDaug=part->daughter(0); jDaug!=0; jDaug=jDaug->next()) {
				jDaug->setP(*decay_.GetDecay(j++));
				double ip(0.);
				ip = getParticleIP(signalpv->getVertex(true),jDaug->getOriginVertex()->getVertex(true),jDaug->getP());
				jDaug->setIP(ip);
				jDaug->smearIP();
				//Now the pileup, we cache the results of the IP smearing first...
				double cachedip = jDaug->getIP();
				double cachedipsmeared = jDaug->getIPSmeared();
				double cachedsigmaip = jDaug->getSigmaIP();
				//The cache for the pileup IP, very dirty but see comment below why...
				double cachedminip = cachedip;
				double cachedminipsmeared = cachedipsmeared;
				double cachedsigmaminip = cachedsigmaip;
				std::vector<RapidVertex>::iterator itrVtx;
				for(itrVtx = pileuppvs_.begin(); itrVtx != pileuppvs_.end(); ++itrVtx) {
					double thisip = getParticleIP((*itrVtx).getVertex(true),jDaug->getOriginVertex()->getVertex(true),jDaug->getP());
					jDaug->setMinIP(thisip);
					jDaug->smearIP();
					if (std::fabs(jDaug->getMinIPSmeared()) < std::fabs(cachedminipsmeared)) {
						cachedminip = jDaug->getMinIP();
						cachedminipsmeared = jDaug->getMinIPSmeared();
						cachedsigmaminip = jDaug->getSigmaMinIP();
					}
				}
				// Use the cached information to set things now... this is not the best coding ever but mandated
				// by the fact that the particle owns the smearing tool...
				jDaug->setIP(cachedip);
				jDaug->setIPSmeared(cachedipsmeared);
				jDaug->setIPSigma(cachedsigmaip);
				jDaug->setMinIP(cachedminip);
				jDaug->setMinIPSmeared(cachedminipsmeared);
				jDaug->setMinIPSigma(cachedsigmaminip);
			}
		}
	}

	return true;
}

double RapidDecay::getParticleIP(ROOT::Math::XYZPoint pv, ROOT::Math::XYZPoint dv, TLorentzVector p) {
	ROOT::Math::XYZVector v1 = pv - dv;
	ROOT::Math::XYZVector lengthv(p.X(), p.Y(), p.Z());
	ROOT::Math::XYZVector v2 = v1 + lengthv;

	ROOT::Math::XYZVector impact = v1.Cross(v2)/sqrt(lengthv.Mag2());

	return sqrt(impact.Mag2());
}

bool RapidDecay::genDecayAccRej() {
	bool passAccRej(true);
	int ntry(0);

	do {
		if(!genDecay(true)) return false;
		passAccRej = runAcceptReject();
		++ntry;

	} while(!passAccRej && ntry<maxgen_);

	if(!passAccRej) {
		if(!suppressAttemptsWarning_) {
			std::cout << "WARNING in RapidDecay::genDecayAccRej : no events found with required kinematics." << std::endl
				  << "                                        this event will not be generated." << std::endl
				  << "                                        further warnings will be suppressed." << std::endl;
			suppressAttemptsWarning_ = true;
		}
		return false;
	}

	return true;
}
