#pragma once


#include <Eigen/core>
#include <Eigen/Dense>
#include "bsp.h"


constexpr double PLANEEPSILON = std::numeric_limits<double>::epsilon() * 100;

struct BspPlane;

struct PlaneTriple
{
    int* p0, * p1, * p2;
};


struct Plane
{
    Eigen::Vector3d Normal;
    double Distance;
};


Eigen::Vector3d IntersectPlanesLU(const Plane& p1, const Plane& p2, const Plane& p3);
Eigen::Vector3d IntersectPlanesLU(const BspPlane& p1_, const BspPlane& p2_, const BspPlane& p3_);
Plane bspplane2plane(const BspPlane& pl);