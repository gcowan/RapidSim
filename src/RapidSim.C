#include <cstdlib>
#include <iostream>
#include <ctime>
#include <iomanip> // For std::fixed and std::setprecision
#include <chrono>  // For std::chrono

#include "TString.h"

#include "RapidAcceptance.h"
#include "RapidConfig.h"
#include "RapidDecay.h"
#include "RapidHistWriter.h"

int rapidSim(const TString mode, const int nEvtToGen, bool saveTree=false, int nToReDecay=0) {

	clock_t t0,t1,t2;

	t0=clock();

	if(!getenv("RAPIDSIM_ROOT")) {
		std::cout << "ERROR in rapidSim : environment variable RAPIDSIM_ROOT is not set" << std::endl
			  << "                    Terminating" << std::endl;
		return 1;
	}

	TString configEnv=getenv("RAPIDSIM_CONFIG");
	if(configEnv!="") {
		std::cout << "INFO in rapidSim : environment variable RAPIDSIM_CONFIG is set" << std::endl
			  << "                   Settings in " << configEnv << " will be used" << std::endl;
	}

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

	if(nToReDecay>0) {
		std::cout << "INFO in rapidSim : re-decay mode is active" << std::endl
			  << "                   Each parent will be re-decayed " << nToReDecay << " times" << std::endl;
	}

	RapidAcceptance* acceptance = config.getAcceptance();

	RapidHistWriter* writer = config.getWriter(saveTree);

	t1=clock();


	// Start the timer
	auto start = std::chrono::high_resolution_clock::now();

	int ngenerated = 0; int nselected = 0;
	for (Int_t n=0; n<nEvtToGen; ++n) {
		writer->setNEvent(n);
		if (!decay->generate()) continue;
		++ngenerated;

		if(acceptance->isSelected()) {
			++nselected;
			writer->fill();
		}

		for (Int_t nrd=0; nrd<nToReDecay; ++nrd) {
			if (!decay->generate(false)) continue;
			++ngenerated;

			if(!acceptance->isSelected()) continue;
			++nselected;

			writer->fill();
		}

		// // Print every 1000 events
		// if (n % 1000 == 0) {
		// 	double percentComplete = (static_cast<double>(n) / nEvtToGen) * 100;
        // 	std::cout << "Event: " << n << " / " << nEvtToGen << " (" << std::fixed << std::setprecision(2) << percentComplete << "% complete)" << std::endl;
		// }

		// Print every 1000 events with percentage completion and time estimation
		if (n % 1000 == 0) {
			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;

			double percentComplete = (static_cast<double>(n) / nEvtToGen) * 100;
			double timePerEvent = elapsed.count() / (n + 1); // time per event
			double estimatedTotalTime = timePerEvent * nEvtToGen; // total estimated time
			double estimatedRemainingTime = estimatedTotalTime - elapsed.count(); // remaining time

			std::cout << "Event: " << n << " / " << nEvtToGen 
					<< " (" << std::fixed << std::setprecision(2) << percentComplete << "% complete)"
					<< ", Time elapsed: " << std::fixed << std::setprecision(2) << elapsed.count() << "s"
					<< ", Estimated remaining time: " << std::fixed << std::setprecision(2) << estimatedRemainingTime << "s"
					<< std::endl;
		}

	}

	writer->save();

	t2=clock();

	std::cout << "INFO in rapidSim : Generated " << ngenerated << std::endl;
	std::cout << "INFO in rapidSim : Selected " << nselected << std::endl;
	std::cout << "INFO in rapidSim : " << (float(t1) - float(t0)) / CLOCKS_PER_SEC << " seconds to initialise." << std::endl;
	std::cout << "INFO in rapidSim : " << (float(t2) - float(t1)) / CLOCKS_PER_SEC << " seconds to generate." << std::endl;

	return 0;
}

int main(int argc, char * argv[])
{
	if (argc < 3) {
		printf("Usage: %s mode numberToGenerate [saveTree=0] [numberToRedecay=0]\n", argv[0]);
		return 1;
	}

	const TString mode = argv[1];
	const int number = static_cast<int>(atof(argv[2]));
	bool saveTree = false;
	int nToReDecay = 0;

	if(argc>3) {
		saveTree = atoi(argv[3]);
	}
	if(argc>4) {
		nToReDecay = atoi(argv[4]);
	}

	int status = rapidSim(mode, number, saveTree, nToReDecay);

	return status;
}
