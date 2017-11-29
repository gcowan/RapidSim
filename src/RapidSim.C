#include <cstdlib>
#include <iostream>

#include "TString.h"

#include "RapidAcceptance.h"
#include "RapidConfig.h"
#include "RapidDecay.h"
#include "RapidHistWriter.h"

int rapidSim(const TString mode, const int nEvtToGen, bool saveTree=false, bool reDecay=false, int nToReDecay=0) {

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

    RapidAcceptance* acceptance = config.getAcceptance();

    RapidHistWriter* writer = config.getWriter(saveTree);

    int ngenerated = 0; int nselected = 0;
    for (Int_t n=0; n<nEvtToGen; ++n) {
        writer->setNEvent(n);     
        if (!decay->generate()) continue;
        ++ngenerated;

        if (reDecay) {
            for (Int_t nrd=0; nrd<nToReDecay; ++nrd) {
                if (!decay->generate(false)) continue;
                ++ngenerated;

                if(!acceptance->isSelected()) continue;
                ++nselected;

                writer->fill();
            }
        } else {
            if(!acceptance->isSelected()) continue;
            ++nselected;

            writer->fill();
        }
    }

    writer->save();

    std::cout << "Generated " << ngenerated << std::endl;
    std::cout << "Selected " << nselected << std::endl;

    return 0;
}

int main(int argc, char * argv[])
{
    if (argc < 3) {
        printf("Usage: %s mode numberToGenerate [saveTree=0] [reDecay=0] [numberToRedecay=0]\n", argv[0]);
        return 1;
    }

    const TString mode = argv[1];
    const int number = atoi(argv[2]);
    bool saveTree = false;
    bool reDecay = false;
    int nToReDecay = 0;

    if(argc>3) {
        saveTree = atoi(argv[3]);
    }
    if(argc>4) {
        reDecay = atoi(argv[4]);
    }
    if(argc>5) {
        nToReDecay = atoi(argv[5]);
    }

    int status = rapidSim(mode, number, saveTree, reDecay, nToReDecay);

    return status;
}
