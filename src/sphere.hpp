#ifndef INC_2DCONTOUR_SPHERE_HPP
#define INC_2DCONTOUR_SPHERE_HPP

#include "pch.h"
#include "shape.hpp"
#include "utils.hpp"

class Sphere : public Shape {
public:
    Sphere(Vector3 c, double _r) : center(c), r(_r) {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = Utils::createIcoSphere(r*0.97, 3, center);
        polyscope::registerSurfaceMesh("Debug"+std::to_string(count++), geometry->vertexPositions, mesh->getFaceVertexList());
    }

    bool isInside(const Vector3& p) noexcept override;
    Vector3 closestSurfacePoint(const Vector3& p) noexcept override;

private:
    Vector3 center;
    double r;
};


#endif //INC_2DCONTOUR_SPHERE_HPP
