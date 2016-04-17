#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

#include "RapidAcceptance.h"
#include "RapidConfig.h"
#include "RapidDecay.h"
#include "RapidHistWriter.h"

int rapidSim(const TString mode, const int nEvtToGen, bool saveTree=false) {

	RapidConfig config;
	if(!config.load(mode)) {
		std::cout << "ERROR in rapidSim : failed to load configuration for decay mode " << mode << std::endl
			  << "                    Terminating" << std::endl;
		return 1;
	}

	RapidDecay* decay = config.getDecay();
	if(!decay) {
		std::cout << "ERROR in rapidSim : failed to setup decay for decay mode " << mode << std::endl
			  << "                    Terminating" << std::endl;
		return 1;
	}

	RapidAcceptance* acceptance = config.getAcceptance();

	RapidHistWriter* writer = config.getWriter(saveTree);

	int ngenerated = 0; int nselected = 0;
	for (Int_t n=0; n<nEvtToGen; ++n) {
		if (!decay->generate()) continue;
		++ngenerated;

		if(!acceptance->inAcceptance()) continue;
		++nselected;

		writer->fill();
	}

	writer->save();

	std::cout << "Generated " << ngenerated << std::endl;
	std::cout << "Selected " << nselected << std::endl;

	return 0;
}

int main(int argc, char * argv[])
{
	if (argc < 3) {
		printf("Usage: %s mode numberToGenerate [saveTree=false]\n", argv[0]);
		return 1;
	}

	const TString mode = argv[1];
	const int number = atoi(argv[2]);
	bool saveTree = false;

	if(argc>3) {
		saveTree = atoi(argv[3]);
	}

	int status = rapidSim(mode, number, saveTree);

	return status;
}
