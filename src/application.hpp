#ifndef COMPLETEMETHOD_APPLICATION_HPP
#define COMPLETEMETHOD_APPLICATION_HPP

#include "pch.h"
#include "triangulate.hpp"
#include "utils.hpp"
#include "meshData.hpp"
#include "GradientDescent.hpp"
#include "ARAP.hpp"
#include "LBFGSopt.hpp"
#include "curvatureProcessing.hpp"
#include "tinyAD.hpp"

enum State {
    CHOOSE_FILE = 0,
    OPTI_1D = 1
};

class Application {
public:

    static void init();
    static void callback0();
    static void callback1();

    static void readOBJ(const std::string& fn);
    static void make2DMesh();
    static void make2DLarge();
    static void show2DMetrics(std::map<size_t, double>& edgeLengths, std::map<size_t, double>& angles);

    // Helper
    static void TraverseDirectory(const std::string& path);
    static void changeState(State newState);

private:
    static void (*callbacks[2])();
    static State state;
    static bool toggle;
    static int highlight;
    static GradientDescent* gd;
    static LBFGSopt* lbfgs;
    static TinyADopt* tinyOpt;
    static ARAP* arapHolder;
    static int index;
    static double x;
    static double y;
    static double z;

    static VertexData<std::array<Vector3,2>>* planeVertexTangent;

    // Mesh
    static std::unique_ptr<ManifoldSurfaceMesh> mesh;
    static std::unique_ptr<VertexPositionGeometry> geometry;

    // 2D mesh
    static std::unique_ptr<ManifoldSurfaceMesh> mesh2D;
    static std::unique_ptr<VertexPositionGeometry> geometry2D;
    static MeshInfo* meshData2D;
    static std::unordered_map<size_t, size_t> matching; // 2D -> 3D

    // Directory
    static std::string currentPath, filename;
    static std::vector<std::string> files, directories;
};


#endif //COMPLETEMETHOD_APPLICATION_HPP
