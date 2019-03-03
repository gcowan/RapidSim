#ifndef RAPIDVERTEX_H
#define RAPIDVERTEX_H

#include <iostream>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

#include "RapidVertexSmear.h"

class RapidVertex {
	public:
		RapidVertex(double x, double y, double z)
			: nPVTracks_(5), vertexTrue_(x, y, z) {smearVertex();}

		ROOT::Math::XYZPoint getVertex(bool truth);

		void setXYZ(double x, double y, double z, RapidVertexSmear* smear = nullptr);

	private:
		void smearVertex();
		void smearVertex(RapidVertexSmear* smear);
		unsigned int nPVTracks_;

		ROOT::Math::XYZPoint vertexTrue_;
		ROOT::Math::XYZPoint vertexSmeared_;
};

#endif
