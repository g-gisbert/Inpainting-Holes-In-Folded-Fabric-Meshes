#ifndef DEF_ARAP
#define DEF_ARAP

#include <iostream>

#include "pch.h"
#include <Eigen/SVD>

using namespace geometrycentral;
using namespace geometrycentral::surface;


class ARAP
{
public:

    ARAP(VertexPositionGeometry& geom);
    ~ARAP();
    void addConstraint(int id, Vector3 value);
    void addNormalConstraint(int id, Vector3 value);
    void solve();

    Eigen::SparseMatrix<double> lapBeltrami;

private:

    void resolveConstraints();
    std::vector<Eigen::Matrix3d> computeRotations();
    void computeNewPositions(std::vector<Eigen::Matrix3d>& R);

    Vector<double> getEigenConstraint(int coord);
    Vector<bool> buildConstraintsVector();
    Eigen::Matrix3d findRotationToFitVectors(Eigen::Vector3d a, Eigen::Vector3d b);

    std::unique_ptr<VertexPositionGeometry> geometryOriginal;
    VertexPositionGeometry& geometryDeformed;

    std::vector<size_t> idNormalConstraints;
    std::vector<Vector3> normalConstraints;
    std::vector<std::pair<size_t, Vector3>> pairedConstraints;
    std::vector<size_t> idConstraints;
    std::array<std::vector<double>, 3> constraints;
    int nConstraints;


    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > *solver;

    int maxIter;
    size_t n; // number of vertices

};

#endif