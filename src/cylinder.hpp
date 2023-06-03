#ifndef INC_2DCONTOUR_CYLINDER_HPP
#define INC_2DCONTOUR_CYLINDER_HPP

#include <utility>

#include "pch.h"
#include "shape.hpp"
#include "utils.hpp"

class Cylinder : public Shape {
public:
    Cylinder(Vector3 c, double _r, double _h, Eigen::Matrix3d  R_) : center(c), r(_r*1.05), h_2(_h/2), R(std::move(R_)) {
        Eigen::Vector3d tmp = R * Eigen::Vector3d{0.0, 1.0, 0.0};
        up = Vector3{tmp[0], tmp[1], tmp[2]};
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = Utils::createCylinder(_r, _h, 32, center, R);
        polyscope::registerSurfaceMesh("DebugCylinder"+std::to_string(count++), geometry->vertexPositions, mesh->getFaceVertexList());
    }

    bool isInside(const Vector3& p) noexcept override;
    Vector3 closestSurfacePoint(const Vector3& p) noexcept override;

private:
    Vector3 center;
    double r;
    double h_2;
    Eigen::Matrix3d R;
    Vector3 up;
};


#endif //INC_2DCONTOUR_CYLINDER_HPP
