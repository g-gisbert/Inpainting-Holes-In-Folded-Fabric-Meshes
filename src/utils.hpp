#ifndef INC_2DCONTOUR_UTILS_HPP
#define INC_2DCONTOUR_UTILS_HPP

#include "pch.h"
#include "ARAP.hpp"
#include "meshData.hpp"
#include "plane.hpp"


struct Triangle {
    Vector3& v1, v2, v3;
    size_t id;

    bool operator==(const Triangle& tri) const {
        return id == tri.id;
    }
};
namespace std {
    template<>
    struct hash<Triangle> {
        std::size_t operator()(const Triangle& tri) const {
            return tri.id;
        }
    };
}

struct Point3d {
    Vector3& v;
    size_t id;

    bool operator==(const Point3d& p) const {
        return id == p.id;
    }
};
namespace std {
    template<>
    struct hash<Point3d> {
        std::size_t operator()(const Point3d& p) const {
            return p.id;
        }
    };
}

struct AABB {
    Vector3 min;
    Vector3 max;
};

namespace Utils {

    double averageBlEdgeLength(const BoundaryLoop& BL, const VertexPositionGeometry& geom);

    //void ARAPinit();
    void ARAP(VertexPositionGeometry& geom, VertexPositionGeometry& geomValues,
                  const MeshInfo& meshData, const std::unordered_map<size_t, size_t>& matching, const double c);

    AABB getAABB(const VertexPositionGeometry& geom);


    std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
        createCylinder(double radius, double height, int N = 32, Vector3 center = Vector3::zero(), Eigen::Matrix3d R = Eigen::Matrix3d::Identity());

    std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
        createIcoSphere(double scale, int level, Vector3 center = Vector3::zero());

    Eigen::Matrix3d eulerAngles(double psi, double phi, double theta);

    Eigen::Vector3d to_eigen(const Vector3& v);
    Vector3 to_geometrycentral(const Eigen::Vector3d& v);
}

#endif //INC_2DCONTOUR_UTILS_HPP
