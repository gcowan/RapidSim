#include "RapidVertexSmearHisto.h"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

RapidVertexSmearHisto::~RapidVertexSmearHisto() {
  while(!histos1D_.empty()) {
    delete histos1D_[histos1D_.size()-1];
    histos1D_.pop_back();
  }

  while(!histos2D_.empty()) {
    delete histos2D_[histos2D_.size()-1];
    histos2D_.pop_back();
  }

  while(!histos3D_.empty()) {
    delete histos3D_[histos3D_.size()-1];
    histos3D_.pop_back();
  }

}


ROOT::Math::XYZPoint RapidVertexSmearHisto::smearVertex(ROOT::Math::XYZPoint vtx, int direction){
  double smearx{0}, smeary{0}, smearz{0};
  unsigned int iHist=0;
  if(histos1D_.size()!=0){
    while(true){//select the right histogram
      if( iHist==thresholds_.size()-1 ) break;
      if( threshVal_ < thresholds_[iHist+1] ) break;
    }
    switch(direction){
    case 0 : smearx = histos1D_.at(iHist)->GetRandom();break;
    case 1 : smeary = histos1D_.at(iHist)->GetRandom();break; 
    case 2 : smeary = histos1D_.at(iHist)->GetRandom();break; 
    default :smearx = histos1D_.at(iHist)->GetRandom();break;
    }
    
  }
  else if(histos2D_.size()!=0){
    while(true){
      if( iHist==thresholds_.size()-1 ) break;
      if( threshVal_ < thresholds_[iHist+1] ) break;
    }
    switch(direction){
    case 0: histos2D_.at(iHist)->GetRandom2(smearx,smeary);break;
    case 1: histos2D_.at(iHist)->GetRandom2(smearx,smearz);break;
    case 2: histos2D_.at(iHist)->GetRandom2(smeary,smearz);break;
    default:histos2D_.at(iHist)->GetRandom2(smearx,smeary);break;
    }
  }
  else if(histos3D_.size()!=0){
    while(true){
      if( iHist==thresholds_.size()-1 ) break;
      if( threshVal_ < thresholds_[iHist+1] ) break;
    }
    histos3D_.at(iHist)->GetRandom3(smearx,smeary,smearz);
  }
  else{
    std::cout<<"Non-configured smearing. Should never come here. BREAK!"<<std::endl;
    assert (0);
  }
  return ROOT::Math::XYZPoint(vtx.X()+smearx,vtx.Y()+smeary,vtx.Z()+smearz);
}


void RapidVertexSmearHisto::init(std::vector<double> thresholds, std::vector<TH1*>histos){
  if(thresholds.size() < histos.size()) {
    std::cout << "WARNING in RapidVertexSmearHisto::init : TH1 : too many histograms provided. Number of histograms should match number of thresholds." << std::endl;
    std::cout << "                                      excess histograms ignored." << std::endl;
    
    while(thresholds.size() < histos.size()) {
      histos.pop_back();
    }
  } else if(thresholds.size() > histos.size()) {
    std::cout << "WARNING in RapidVertexSmearHisto::init : TH1 : too few histograms provided. Number of histograms should match number of thresholds." << std::endl;
    std::cout << "                                      excess thresholds ignored." << std::endl;
    
    while(thresholds.size() > histos.size()) {
      thresholds.pop_back();
    }
  }
  
  thresholds_ = thresholds;
  histos1D_ = histos;
}
void RapidVertexSmearHisto::init(std::vector<double> thresholds, std::vector<TH2*>histos){
  if(thresholds.size() < histos.size()) {
    std::cout << "WARNING in RapidVertexSmearHisto::init : TH2 : too many histograms provided. Number of histograms should match number of thresholds." << std::endl;
    std::cout << "                                      excess histograms ignored." << std::endl;
    
    while(thresholds.size() < histos.size()) {
      histos.pop_back();
    }
  } else if(thresholds.size() > histos.size()) {
    std::cout << "WARNING in RapidVertexSmearHisto::init : TH2 : too few histograms provided. Number of histograms should match number of thresholds." << std::endl;
    std::cout << "                                      excess thresholds ignored." << std::endl;
    
    while(thresholds.size() > histos.size()) {
      thresholds.pop_back();
    }
  }
  
  thresholds_ = thresholds;
	histos2D_ = histos;
}
void RapidVertexSmearHisto::init(std::vector<double> thresholds, std::vector<TH3*>histos){
    if(thresholds.size() < histos.size()) {
      std::cout << "WARNING in RapidVertexSmearHisto::init : TH3 : too many histograms provided. Number of histograms should match number of thresholds." << std::endl;
      std::cout << "                                      excess histograms ignored." << std::endl;
      
      while(thresholds.size() < histos.size()) {
			histos.pop_back();
      }
    } else if(thresholds.size() > histos.size()) {
      std::cout << "WARNING in RapidVertexSmearHisto::init : TH3 : too few histograms provided. Number of histograms should match number of thresholds." << std::endl;
      std::cout << "                                      excess thresholds ignored." << std::endl;
      
      while(thresholds.size() > histos.size()) {
	thresholds.pop_back();
      }
    }
    
    thresholds_ = thresholds;
    histos3D_ = histos;
}
