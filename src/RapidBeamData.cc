#include "RapidBeamData.h"

#include <cstdlib>
#include <fstream>
#include <iostream>

RapidBeamData* RapidBeamData::instance_=0;

RapidBeamData* RapidBeamData::getInstance() {
	if(!instance_) {
		instance_ = new RapidBeamData();
		TString path;
		path=getenv("RAPIDSIM_CONFIG");
		if(path!="") instance_->loadData(path+"/config/beam.dat");

		path=getenv("RAPIDSIM_ROOT");
		instance_->loadData(path+"/config/beam.dat");
	}
	return instance_;
}


void RapidBeamData::loadData(TString file) {
	std::cout << "INFO in RapidBeamData::loadData : loading beam data from " << file << std::endl;

	std::ifstream fin;
	fin.open(file, std::ifstream::in);
	TString buffer;
	buffer.ReadLine(fin);//ignore title line

	while(fin.good()) {
		buffer.ReadToken(fin);
		pileup_ = buffer.Atoi();
		buffer.ReadToken(fin);
		sigmaxy_ = buffer.Atof();
		buffer.ReadToken(fin);
		sigmaz_ = buffer.Atof();

		if(!fin.good()) break;

	}
	fin.close();
}


