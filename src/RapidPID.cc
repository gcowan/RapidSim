#include "RapidPID.h"

#include <iostream>

RapidPID::~RapidPID() {
	std::map<unsigned int, TH3D*>::iterator itr = pidHists_.begin();
	while (itr != pidHists_.end()) {
		delete itr->second;
		pidHists_.erase(itr++);
	}

	std::map<unsigned int, std::map<unsigned int, TH1D*>*>::iterator itr2 = cachedPIDHists_.begin();
	while (itr2 != cachedPIDHists_.end()) {
		std::map<unsigned int, TH1D*>::iterator itr3 = (*itr2).second->begin();
		while (itr3 != (*itr2).second->end()) {
			delete itr3->second;
			(*itr2).second->erase(itr3++);
		}
		delete itr2->second;
		cachedPIDHists_.erase(itr2++);
	}
}

double RapidPID::getPID(unsigned int id, double p, double eta) {
	//first check that we have a histogram for the given particle ID
	if(pidHists_.find(id)==pidHists_.end()) {
		std::cout << "WARNING in RapidPID::getPID : PID histogram not set for " << name_ << " " << id << std::endl;
		std::cout << "                              returning 0" << std::endl;
		return 0.;
	}
	TH3D* h = pidHists_.at(id);
	if(!h) {
		std::cout << "WARNING in RapidPID::getPID : PID histogram not set for " << name_ << " " << id << std::endl;
		std::cout << "                              returning 0" << std::endl;
		return 0.;
	}

	//if p and/or eta is out of range then set to the limit
	limitRangePEta(id,p,eta);

	//now get the bin (X,Y) and check if we've already cached the projected TH1
	unsigned int nX   = h->GetNbinsX()+2;
	unsigned int nY   = h->GetNbinsY()+2;

	unsigned int bin  = h->FindBin(p,eta,-1.) % (nX*nY);
	unsigned int binX = bin%nX;
	unsigned int binY = ((bin-binX)/nX)%nY;

	if(cachedPIDHists_.find(id)==cachedPIDHists_.end()) {
		cachedPIDHists_.insert(std::pair<unsigned int, std::map<unsigned int, TH1D*>*>(id,new std::map<unsigned int, TH1D*>()));
	}

	std::map<unsigned int, TH1D*>* cachedBins = cachedPIDHists_.at(id);

	//if bin isn't cached let's project it and cache it now
	if(cachedBins->find(bin)==cachedBins->end()) {
		TString hname = "cachedPID"; hname+=name_; hname+="_"; hname+=id; hname+="_"; hname+=bin;
		cachedBins->insert(std::pair<unsigned int,TH1D*>(bin,pidHists_.at(id)->ProjectionZ(hname, binX, binX+1, binY, binY+1)));
	}

	TH1D* hist = cachedBins->at(bin);

	return hist->GetRandom();
}

void RapidPID::addPID(unsigned int id, TH3D* hist) {
	pidHists_.insert(std::pair<unsigned int,TH3D*>(id,hist));
	setupRangePEta(id);
}

void RapidPID::setupRangePEta(unsigned int id) {
	TH3D* hist = pidHists_.at(id);
	if(!hist) return;

	TH1D* hist_x = hist->ProjectionX();
	TH1D* hist_y = hist->ProjectionY();
	hist_x->Sumw2(false);
	hist_y->Sumw2(false);

	maxP_[id]   = hist_x->FindLastBinAbove (0);
	minEta_[id] = hist_y->FindFirstBinAbove(0);
	maxEta_[id] = hist_y->FindLastBinAbove (0);
}

void RapidPID::limitRangePEta(unsigned int id, double& p, double& eta) {
	if(p > maxP_[id]) p = maxP_[id];
	if(eta < minEta_[id]) eta = minEta_[id];
	if(eta > maxEta_[id]) eta = maxEta_[id];
}
