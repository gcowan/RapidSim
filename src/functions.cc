#include "functions.h"

bool generateEvent(TLorentzVector& head, TGenPhaseSpace& event, double* masses , int np, TRandom& ran, int m_maxgen ){

    /* TLorentzVector head Particle to decay
       TGenPhaseSpace decay generator
       double* masses array of output particles
       int m_maxgen number to try  
       */

    // check decay kinematics valid
    bool isok = event.SetDecay(head, np, masses);
    if (!isok) return false;

    // make an event
    int ntoGen = 0; bool accept = false;
    while (ntoGen < m_maxgen && accept == false){
        Double_t weight = event.Generate();
        accept = weight > ran.Uniform();
        ++ntoGen;
    } // while

    return accept;

}


double pick(RooDataSet * data ,TRandom& ran, std::string var_name){
    int entry = int(data->numEntries() * ran.Uniform());
    const RooArgSet* row = data->get(entry); 
    double value = row->getRealValue(var_name.c_str()); 
    return value;
}

std::string varname(std::string header, std::string var){
    return (header + "_" + var);
}

RooGounarisSakurai* rooGS(RooRealVar& m, double mean, double gamma,
        double thespin, double m1,double m2, std::string name){ 
    // rho+ -> pi0 pi+ decay
    RooRealVar* m0 = new RooRealVar(varname(name,std::string("m0")).c_str(),varname(name,std::string("m0")).c_str(), mean);
    RooRealVar* g0 = new RooRealVar(varname(name,std::string("g0")).c_str(),varname(name,std::string("g0")).c_str(),gamma);
    RooRealVar* spin = new RooRealVar(varname(name,std::string("spin")).c_str(),varname(name,std::string("spin")).c_str(),thespin);
    RooRealVar* radius = new RooRealVar(varname(name,std::string("radius")).c_str(),varname(name,std::string("radius")).c_str(), barrierFactor); // not used
    RooRealVar* ma = new RooRealVar(varname(name,std::string("ma")).c_str(),varname(name,std::string("ma")).c_str(),m1);
    RooRealVar* mb = new RooRealVar(varname(name,std::string("mb")).c_str(),varname(name,std::string("mb")).c_str(),m2);   

    return new RooGounarisSakurai(name.c_str(),name.c_str(), m,*m0,*g0,*spin, *radius,*ma,*mb);
}


RooRelBreitWigner* rooBW(RooRealVar& m,  double mean, double gamma,
        double thespin, double m1,double m2, std::string name){

    RooRealVar* m0 = new RooRealVar(varname(name,std::string("m0")).c_str(),varname(name,std::string("m0")).c_str(), mean);
    RooRealVar* g0 = new RooRealVar(varname(name,std::string("g0")).c_str(),varname(name,std::string("g0")).c_str(),gamma);
    RooRealVar* spin= new RooRealVar(varname(name,std::string("spin")).c_str(),varname(name,std::string("spin")).c_str(),thespin);
    RooRealVar* radius= new RooRealVar(varname(name,std::string("radius")).c_str(),varname(name,std::string("radius")).c_str(),barrierFactor); // not used
    RooRealVar* ma= new RooRealVar(varname(name,std::string("ma")).c_str(),varname(name,std::string("ma")).c_str(), m1);
    RooRealVar* mb= new RooRealVar(varname(name,std::string("mb")).c_str(),varname(name,std::string("mb")).c_str(), m2);   

    return new RooRelBreitWigner(name.c_str(),name.c_str(), m,*m0,*g0,*radius,*ma,*mb,*spin);
}


RooGounarisSakurai* createRhoPlus(RooRealVar& m, std::string name){
    return rooGS(m,0.77511,0.1491,1, 0.13957, 0.1349766 , name);
}

RooRelBreitWigner* createPhiMassPdf(RooRealVar& m, std::string name){
    return rooBW(m,mphi ,0.001491,1,mK ,mK ,name ); 
}

RooRelBreitWigner* createKstarMassPdf(RooRealVar& m, std::string name){
    return rooBW(m,0.89166 ,0.0508,1,mK ,mpi ,name ); 
}


RooRelBreitWigner* createChi0MassPdf(RooRealVar& m, std::string name){
    return rooBW(m,mchic0 ,10.5e-3,1,mpi,mpi,name ); 
}

RooRelBreitWigner* createChi1MassPdf(RooRealVar& m, std::string name){
    return rooBW(m,mchic ,0.31e-3,1,mpi,mpi,name ); 
}

RooRelBreitWigner* createChi2MassPdf(RooRealVar& m, std::string name){
    return rooBW(m,mchic2 ,1.93e-3,1,mpi,mpi,name ); 
}

RooRelBreitWigner* createpsi2MassPdf(RooRealVar& m, std::string name){
    return rooBW(m,mPsi ,0.3e-3,1,mpi,mpi,name ); 
}


TLorentzVector genB(TRandom ran, TH1F* ptHisto, TH1F* etaHisto, double m){
    TLorentzVector vec;
    double phi = ran.Uniform(0,2*TMath::Pi());
    vec.SetPtEtaPhiM(ptHisto->GetRandom(),etaHisto->GetRandom(),phi,m);
    return vec;
}

double resSlope(double p) {
    // std::cout << p << std::endl;
    return(sqrt(pow(6.2e-5,2) + pow(2.1e-3/p,2)));  
}

TLorentzVector toFourVector(const TVector3& vec, double m) {
    return TLorentzVector(vec, TMath::Sqrt(m*m + vec.Mag2()));
}

TLorentzVector reassignMass(const TLorentzVector& vec, double mass){
    const TVector3 threevec = vec.Vect();
    return toFourVector(threevec, mass);
}

TVector3 smearedVec(double plus_x, double plus_y, double plus_z, TGraphErrors* dGraph, TRandom& ran) {
    double kp, kptx, kpty, norm, smear;
    kp = sqrt(plus_x*plus_x + plus_y*plus_y + plus_z*plus_z );
    kptx = plus_x/plus_z;
    kpty = plus_y/plus_z;
    norm = sqrt(1 + kptx*kptx + kpty*kpty);
    smear = 1.0*ran.Gaus(0,1)*dGraph->Eval(1000*kp)*kp;
    kp += smear;

    // smear the slopes
    double slope_smear = resSlope(kp); 
    kptx += slope_smear*ran.Gaus(1,0);
    kpty += slope_smear*ran.Gaus(1,0);
    norm = sqrt(1 + kptx*kptx + kpty*kpty);

    plus_x = kptx*kp/norm;
    plus_y = kpty*kp/norm;
    plus_z = kp/norm;
    return TVector3(plus_x,plus_y, plus_z);

}

TLorentzVector smearedVec(TLorentzVector& vec, TGraphErrors* dGraph, TRandom& ran){

    TVector3 threeVec = smearedVec(vec.Px(), vec.Py(), vec.Pz(),dGraph,ran);
    return toFourVector(threeVec, vec.M());
}

bool inAcceptance(TLorentzVector& vec){

    if (TMath::Abs(vec.Px()/vec.Pz()) > 0.3) return false;
    if (TMath::Abs(vec.Py()/vec.Pz()) > 0.25) return false;  
    if (sqrt(pow(vec.Px()/vec.Pz(),2) + pow(vec.Py()/vec.Pz(),2)) <0.01) return false;

    return true;
}

TLorentzVector magnetKick(TLorentzVector& vec, double charge){

    TVector3 threeVec  =  vec.Vect();

    // extrapolate to magnet centre  

    // kick 
    double p = vec.P();
    double px =  vec.Px() + ptkick*charge ;
    double py = vec.Py();
    double pz = sqrt(p*p - px*px - py*py);

    TVector3 threevec = TVector3(px,py,pz) ;

    return toFourVector(threevec, vec.M()); 

}

bool inDownstream(TLorentzVector& vec, int charge){

    TLorentzVector newvec = magnetKick(vec,charge);

    // position at magnet centre
    double xMag = zC*vec.Px()/vec.Pz(); 
    double yMag = zC*vec.Py()/vec.Pz();  

    double xTracker = xMag +  (newvec.Px()*(zTracker - zC)/newvec.Pz()); 
    double yTracker = yMag +  (newvec.Py()*(zTracker - zC)/newvec.Pz());   

    //if (TMath::Abs(xTracker) > xSizeTracker) std::cout << "out of tracker " << std::endl;

    return (TMath::Abs(xTracker) < xSizeTracker && TMath::Abs(xTracker) > xMinTracker 
            &&  TMath::Abs(yTracker) < ySizeTracker) && TMath::Abs(yTracker) > yMinTracker ;  
}

int pdgCode(TString part) {
	//leptons
	if(part=="e-")         return  11;
	if(part=="e+")         return -11;
	if(part=="nue")        return  12;
	if(part=="anti-nue")   return -12;
	if(part=="mu-")        return  13;
	if(part=="mu+")        return -13;
	if(part=="numu")       return  14;
	if(part=="anti-numu")  return -14;
	if(part=="tau-")       return  15;
	if(part=="tau+")       return -15;
	if(part=="nutau")      return  16;
	if(part=="anti-nutau") return -16;

	//light mesons
	if(part=="pi0")       return      111;
	if(part=="pi+")       return      211;
	if(part=="pi-")       return     -211;
	if(part=="eta")       return      221;
	if(part=="rho0")      return      113;
	if(part=="rho+")      return      213;
	if(part=="rho-")      return     -213;
	if(part=="omega")     return      223;
	
	//strange mesons
	if(part=="KL")        return      130;
	if(part=="KS")        return      310;
	if(part=="K0")        return      311;
	if(part=="K0b")       return     -311;
	if(part=="K+")        return      321;
	if(part=="K-")        return     -321;
	if(part=="K*0")       return      313;
	if(part=="K*0b")      return     -313;
	if(part=="K*+")       return      323;
	if(part=="K*-")       return     -323;

	//ssbar
	if(part=="eta'")      return      331;
	if(part=="phi")       return      333;

	//charm mesons
	if(part=="D+")        return      411;
	if(part=="D-")        return     -411;
	if(part=="D0")        return      421;
	if(part=="D0b")       return     -421;
	if(part=="Ds+")       return      431;
	if(part=="Ds-")       return     -431;
	if(part=="D*+")       return      413;
	if(part=="D*-")       return     -413;
	if(part=="D*0")       return      423;
	if(part=="D*0b")      return     -423;
	if(part=="Ds*+")      return      433;
	if(part=="Ds*-")      return     -433;

	//ccbar
	if(part=="etac")      return      441;
	if(part=="chic0")     return    10441;
	if(part=="Jpsi")      return      443;
	if(part=="chic1")     return    20443;
	if(part=="psi2S")     return   100443;
	if(part=="chic2")     return      445;
	
	//bottom mesons
	if(part=="B0")        return      511;
	if(part=="B0b")       return     -511;
	if(part=="B+")        return      521;
	if(part=="B-")        return     -521;
	if(part=="Bs0")       return      531;
	if(part=="Bs0b")      return     -531;
	if(part=="Bc+")       return      541;
	if(part=="Bc-")       return     -541;
	if(part=="B*0")       return      513;
	if(part=="B*0b")      return     -513;
	if(part=="B*+")       return      523;
	if(part=="B*-")       return     -523;
	if(part=="Bs*0")      return      533;
	if(part=="Bs*0b")     return     -533;

	//bbbar
	if(part=="etab")      return      551;
	if(part=="Upsilon1S") return      553;
	if(part=="Upsilon2S") return   100553;
	if(part=="Upsilon3S") return   200553;
	if(part=="Upsilon4S") return   300553;

	//baryons
	if(part=="n")         return     2112;
	if(part=="p")         return     2212;
	if(part=="Lambda")    return     3122;
	if(part=="Sigma-")    return     3112;
	if(part=="Sigma0")    return     3212;
	if(part=="Sigma+")    return     3222;
	if(part=="Lambda_c")  return     4122;
	if(part=="Sigma_c0")  return     4112;
	if(part=="Sigma_c+")  return     4212;
	if(part=="Sigma_c++") return     4222;
	if(part=="Lambda_b")  return     5122;
	if(part=="Sigma_b-")  return     5112;
	if(part=="Sigma_b+")  return     5222;

	if(part=="anti-n")         return    -2112;
	if(part=="anti-p")         return    -2212;
	if(part=="anti-Lambda")    return    -3122;
	if(part=="anti-Sigma+")    return    -3112;
	if(part=="anti-Sigma0")    return    -3212;
	if(part=="anti-Sigma-")    return    -3222;
	if(part=="anti-Lambda_c")  return    -4122;
	if(part=="anti-Sigma_c0")  return    -4112;
	if(part=="anti-Sigma_c-")  return    -4212;
	if(part=="anti-Sigma_c--") return    -4222;
	if(part=="anti-Lambda_b")  return    -5122;
	if(part=="anti-Sigma_b+")  return    -5112;
	if(part=="anti-Sigma_b-")  return    -5222;

	std::cout << "WARNING in pdgCode : unknown particle " << part << " will be treated as massless" << std::endl;
	return 0;
}

double getMass(int pdgCode) {
	switch(pdgCode) {
		case       11:
		case      -11:
			return me;
		case       13:
		case      -13:
			return mmu;
		case       15:
		case      -15:
			return mtau;
		case      111:
			return mpi0;
		case      211:
		case     -211:
			return mpi;
		case      221:
			return meta;
		case      113:
			return mrho0;
		case      213:
		case     -213:
			return mrho;
		case      223:
			return momega;
		case      130:
		case      310:
		case      311:
		case     -311:
			return mK0;
		case      321:
		case     -321:
			return mK;
		case      313:
		case     -313:
			return mKstar;
		case      323:
		case     -323:
			return mKstar0;
		case      331:
			return metapr;
		case      333:
			return mphi;
		case      411:
		case     -411:
			return mDplus;
		case      421:
		case     -421:
			return mD0;
		case      431:
		case     -431:
			return mDs;
		case      413:
		case     -413:
			return mDstar;
		case      423:
		case     -423:
			return mD0star;
		case      433:
		case     -433:
			return mDsstar;
		case      441:
			return metac;
		case    10441:
			return mchic0;
		case      443:
			return mJpsi;
		case    20443:
			return mchic;
		case   100443:
			return mPsi;
		case      445:
			return mchic2;
		case      511:
		case     -511:
			return mBd;
		case      521:
		case     -521:
			return mB;
		case      531:
		case     -531:
			return mBs;
		case      541:
		case     -541:
			return mBc;
		case      513:
		case     -513:
			return mB0star;
		case      523:
		case     -523:
			return mBstar;
		case      533:
		case     -533:
			return mBsstar;
		case      553:
			return mUp1S;
		case      551:
			return metab;
		case   100553:
			return mUp2S;
		case   200553:
			return mUp3S;
		case   300553:
			return mUp4S;
		case     2112:
		case    -2112:
			return mn;
		case     2212:
		case    -2212:
			return mp;
		case     3122:
		case    -3122:
			return mL;
		case     3112:
		case    -3112:
			return mSm;
		case     3212:
		case    -3212:
			return mS0;
		case     3222:
		case    -3222:
			return mSp;
		case     4122:
		case    -4122:
			return mLc;
		case     4112:
		case    -4112:
			return mSc0;
		case     4212:
		case    -4212:
			return mScp;
		case     4222:
		case    -4222:
			return mScpp;
		case     5122:
		case    -5122:
			return mLb;
		case     5112:
		case    -5112:
			return mSbm;
		case     5222:
		case    -5222:
			return mSbp;
	}
	return 0.;
}
