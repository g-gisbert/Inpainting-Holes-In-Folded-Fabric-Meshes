
#ifndef DEF_CURVATUREPROCESSING
#define DEF_CURVATUREPROCESSING
 
#include "pch.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

struct curvatureData {

    float k1;
    Vector3 u1;
    float k2;
    Vector3 u2;
    Vector3 n;
    int pointType;
    int boundary;
    float debug;
    Vector3 dDebug;
    float H;
    float dHx;
    float dHy;
    float dHz;
    float dk1x;
    float dk1y;
    float dk1z;
    float dk2x;
    float dk2y;
    float dk2z;
    EdgeData<float> edgeDebug;
    EdgeData<Vector3> edgeDebugDerivative;
    
    curvatureData(){};
};
 
class curvatureProcessing
{
    public:
 
    curvatureProcessing();

    static VertexData<curvatureData> calculateCurvature(VertexPositionGeometry& geometry, float r_coeff);

};
 
#endif