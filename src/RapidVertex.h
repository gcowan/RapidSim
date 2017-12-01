#ifndef RAPIDVERTEX_H
#define RAPIDVERTEX_H

#include <iostream>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

class RapidVertex {
	public:
		RapidVertex(double x, double y, double z)
			: vertexTrue_(x, y, z) {smearVertex();}

        std::pair<ROOT::Math::XYZPoint, ROOT::Math::XYZPoint> getVertex() {return vertex_;};

	private:
        void smearVertex();
        
        ROOT::Math::XYZPoint vertexTrue_;
        ROOT::Math::XYZPoint vertexSmeared_;
        std::pair<ROOT::Math::XYZPoint, ROOT::Math::XYZPoint> vertex_;
};

#endif
