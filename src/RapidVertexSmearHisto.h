#ifndef RAPIDVERTEXSMEARHISTO_H
#define RAPIDVERTEXSMEARHISTO_H 1

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "RapidVertexSmear.h"


class RapidVertexSmearHisto : public RapidVertexSmear {
        ///Simple class to smear Vertex position by a TH1, TH2 or TH3.
        //variable 'direction' controls which direction to smear
        // TH1: 0 == x,  1== y,  2== z
        // TH2: 0 == xy, 2== xz, 3== yz
        // TH3: only 0==xyz
        //thresholds allow for possible extra variable dependence, e.g. decay time.
        
        public:
                RapidVertexSmearHisto(std::vector<double> thresholds, std::vector<TH1*> histos){ init(thresholds, histos); };
	        RapidVertexSmearHisto(std::vector<double> thresholds, std::vector<TH2*> histos){ init(thresholds, histos); };
	        RapidVertexSmearHisto(std::vector<double> thresholds, std::vector<TH3*> histos){ init(thresholds, histos); };
	        ~RapidVertexSmearHisto();

		ROOT::Math::XYZPoint smearVertex(ROOT::Math::XYZPoint vtx, int direction);
		void SetThresholdVariableValue(double val){threshVal_=val;}
        private:
		void init(std::vector<double> thresholds, std::vector<TH1*>histos);
		void init(std::vector<double> thresholds, std::vector<TH2*>histos);
		void init(std::vector<double> thresholds, std::vector<TH3*>histos);
		std::vector<TH1*> histos1D_;
		std::vector<TH2*> histos2D_;
		std::vector<TH3*> histos3D_;
		std::vector<double>thresholds_;
		double threshVal_;
};



#endif
