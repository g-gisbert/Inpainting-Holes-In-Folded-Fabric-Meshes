#ifndef INC_2DCONTOUR_GRADIENTDESCENT_HPP
#define INC_2DCONTOUR_GRADIENTDESCENT_HPP

#include "pch.h"
#include "meshData.hpp"
#include "shape.hpp"
#include "sphere.hpp"
#include "cylinder.hpp"
#include "utils.hpp"

class GradientDescent {
public:
    explicit GradientDescent(MeshInfo& _meshData) : alpha(0.1), beta(0.0), gamma(0.001), lambda(1.0), mu(0.01), timeStep(0.5),
        meshData(_meshData),
        grad(std::vector<Vector3>(meshData.mesh.nVertices(), Vector3::zero())), iter(0) {
        //objects.push_back(std::make_unique<Sphere>(Vector3{0.0, 0.5, 0.0}, 0.509));
        //objects.push_back(std::make_unique<Sphere>(Vector3{-0.660382, 0.749638, 0.0}, 1.0));
        //objects.push_back(std::make_unique<Cylinder>(Vector3{-3.65421, 0.652517, 0.0}, 1.0, 6.0, Utils::eulerAngles(0.0, -88.1/180.0*M_PI, 0.0)));
        //objects.push_back(std::make_unique<Cylinder>(Vector3{0.114316, 3.64177, 0.0}, 1.0, 6.0, Utils::eulerAngles(0, -15.0/180.0*M_PI, 0)));

        AABB aabb = Utils::getAABB(meshData.geom);
        gridSize = max(max(aabb.max.y - aabb.min.y, aabb.max.z - aabb.min.z), abs(aabb.max.x - aabb.min.x)) / 10.0;
    }

    void step() noexcept;
    void update() noexcept;
    void updateSpatialHash() noexcept;
    [[nodiscard]] int hash(const Vector3& v) const noexcept;
    void detectCollisions() noexcept;

    float alpha;
    float beta;
    float gamma;
    float lambda;
    float mu;
    float timeStep;

private:
    std::vector<Vector3> previousState;
    MeshInfo& meshData;
    std::vector<Vector3> grad;
    std::vector<Vertex> collidingVertices;
    std::vector<Vector3> anchorPoints;
    std::vector<std::unique_ptr<Shape>> objects;

    //std::unordered_map<int, std::unordered_set<Triangle>> spatialHashing;
    std::unordered_map<int, std::vector<Point3d>> spatialHashing;
    double gridSize;
    int iter;
};


#endif //INC_2DCONTOUR_GRADIENTDESCENT_HPP
