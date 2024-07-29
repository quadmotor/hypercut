
#include "util.h"





// double check
Eigen::Vector3d IntersectPlanesLU(const Plane& p1, const Plane& p2, const Plane& p3)
{
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    A << p1.Normal.x(), p1.Normal.y(), p1.Normal.z(),
        p2.Normal.x(), p2.Normal.y(), p2.Normal.z(),
        p3.Normal.x(), p3.Normal.y(), p3.Normal.z();

    b << p1.Distance, p2.Distance, p3.Distance;

    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
    return x;
}


Plane bspplane2plane(const BspPlane& pl)
{
    Plane rawpl;
    rawpl.Normal << pl.a, pl.b, pl.c;
    const double len = rawpl.Normal.norm();
    rawpl.Distance = -pl.d / len;
    rawpl.Normal /= len;

    return rawpl;
}


Eigen::Vector3d IntersectPlanesLU(const BspPlane& p1_, const BspPlane& p2_, const BspPlane& p3_)
{
    auto p1 = bspplane2plane(p1_);
    auto p2 = bspplane2plane(p2_);
    auto p3 = bspplane2plane(p3_);

    return IntersectPlanesLU(p1, p2, p3);
}
