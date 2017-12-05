#ifndef RAPIDEXTERNALEVTGEN_H
#define RAPIDEXTERNALEVTGEN_H

#include "TString.h"

#include "RapidExternalGenerator.h"

#ifdef RAPID_EVTGEN
class EvtGen;
#endif

class RapidExternalEvtGen : public RapidExternalGenerator {
	public:
		RapidExternalEvtGen()
#ifdef RAPID_EVTGEN
			: evtGen_(0)
#else
			: suppressWarning_(false)
#endif
			{}
		virtual bool decay(std::vector<RapidParticle*>& parts);
		virtual bool setup();
		
		bool setupGenerator();
		void writeDecFile(TString fname, std::vector<RapidParticle*>& parts);

		static TString getEvtGenName(int id);
		static TString getEvtGenConjName(int id);

	private:
#ifdef RAPID_EVTGEN
		EvtGen* evtGen_;
#else
		bool suppressWarning_;
#endif
		TString decFileName_;
};

#endif

