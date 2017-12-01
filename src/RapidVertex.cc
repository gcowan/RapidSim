#include "RapidVertex.h"

#include <iostream>

#include "TMath.h"

void RapidVertex::smearVertex() {
    double xS = 0.010817 + 0.03784*TMath::Exp(-0.0815*vertexTrue_.X());
    double yS = 0.010817 + 0.03784*TMath::Exp(-0.0815*vertexTrue_.Y());
    double zS = 0.010817 + 0.03784*TMath::Exp(-0.0815*vertexTrue_.Z());
    vertexSmeared_ = ROOT::Math::XYZPoint(xS, yS, zS);
    vertex_ = std::make_pair(vertexTrue_, vertexSmeared_);
}

