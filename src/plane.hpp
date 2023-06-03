#ifndef DEF_PLANE
#define DEF_PLANE
 
#include "pch.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

 
struct Plane
{

    Plane(Vector3 _n);
    Plane(Vector3 _n, double _shift);
    Vector3 projection(Vector3 p) const;
    
    Vector3 n;
    double shift;

};

 
#endif