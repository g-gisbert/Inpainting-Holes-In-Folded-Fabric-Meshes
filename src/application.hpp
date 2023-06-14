#ifndef COMPLETEMETHOD_APPLICATION_HPP
#define COMPLETEMETHOD_APPLICATION_HPP

#include "pch.h"
#include "triangulate.hpp"
#include "utils.hpp"
#include "meshData.hpp"
#include "GradientDescent.hpp"
#include "ARAP.hpp"

class Application {
public:

    static void init(const std::string& df);

    static void readOBJ(const std::string& fn);
    static void make2DMesh(const std::string& fn);
    static void show2DMetrics(std::map<size_t, double>& edgeLengths, std::map<size_t, double>& angles);

    static void callback();

private:
    static bool toggle;
    static bool bendingEnergy;
    static bool radioOption;
    static bool enableCollisions;
    static GradientDescent* gd;
    static int index;


    // Mesh
    static std::unique_ptr<ManifoldSurfaceMesh> mesh;
    static std::unique_ptr<VertexPositionGeometry> geometry;

    // 2D mesh
    static std::unique_ptr<ManifoldSurfaceMesh> mesh2D;
    static std::unique_ptr<VertexPositionGeometry> geometry2D;
    static MeshInfo* meshData2D;
    static std::unordered_map<size_t, size_t> matching; // 2D -> 3D

};


#endif //COMPLETEMETHOD_APPLICATION_HPP
