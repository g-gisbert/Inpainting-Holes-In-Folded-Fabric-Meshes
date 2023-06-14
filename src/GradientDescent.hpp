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
    explicit GradientDescent(MeshInfo& _meshData) : alpha(0.1), beta(0.02), mu(0.01), timeStep(0.5), nObjects(0),
        meshData(_meshData),
        grad(std::vector<Vector3>(meshData.mesh.nVertices(), Vector3::zero())) {
        parseObjects();
        AABB aabb = Utils::getAABB(meshData.geom);
        gridSize = max(max(abs(aabb.max.y - aabb.min.y), abs(aabb.max.z - aabb.min.z)), abs(aabb.max.x - aabb.min.x)) / 10.0;
    }

    void step(const bool bendingEnergy, const bool enableCollisions) noexcept;
    void update() noexcept;
    void updateSpatialHash() noexcept;
    [[nodiscard]] int hash(const Vector3& v) const noexcept;
    void detectCollisions() noexcept;

    float alpha; // edge lengths
    float beta; // bending
    float mu; // collisions
    float timeStep;

    int nObjects;

private:
    void parseObjects();

    std::vector<Vector3> previousState;
    MeshInfo& meshData;
    std::vector<Vector3> grad;
    std::vector<Vertex> collidingVertices;
    std::vector<Vector3> anchorPoints;
    std::vector<std::unique_ptr<Shape>> objects;

    std::unordered_map<int, std::vector<Point3d>> spatialHashing;
    double gridSize;
};


#endif //INC_2DCONTOUR_GRADIENTDESCENT_HPP
