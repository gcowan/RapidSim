#ifndef RAPIDEXTERNALGENERATOR_H
#define RAPIDEXTERNALGENERATOR_H

#include <vector>

#include "RapidParticle.h"

class RapidExternalGenerator {
	public:
		virtual ~RapidExternalGenerator() {}

		virtual bool decay(std::vector<RapidParticle*>& /*parts*/)=0;
		virtual bool setup()=0;
};

#endif

