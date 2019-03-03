#ifndef RAPIDVERTEXSMEAR_H
#define RAPIDVERTEXSMEAR_H 1

#include "Math/Point3D.h"

class RapidVertexSmear {
        public:
                virtual ~RapidVertexSmear() {}

		virtual ROOT::Math::XYZPoint smearVertex(ROOT::Math::XYZPoint vtx, int direction)=0;
};


#endif
