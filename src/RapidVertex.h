#ifndef RAPIDVERTEX_H
#define RAPIDVERTEX_H

#include <iostream>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

class RapidVertex {
	public:
		RapidVertex(double x, double y, double z)
			: ntracks_(4), vertexTrue_(x, y, z) {smearVertex();}

		ROOT::Math::XYZPoint getVertex(bool truth);

		void setXYZ(double x, double y, double z);
		void setNtracks(unsigned int ntracks) { ntracks_ = ntracks; smearVertex();}

	private:
		void smearVertex();

		unsigned int ntracks_;
		ROOT::Math::XYZPoint vertexTrue_;
		ROOT::Math::XYZPoint vertexSmeared_;
};

#endif
