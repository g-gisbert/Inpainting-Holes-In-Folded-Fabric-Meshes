#include "ARAP.hpp"


#include <functional>
#include <vector>


using namespace geometrycentral;
using namespace geometrycentral::surface;


ARAP::ARAP(VertexPositionGeometry& geom) : geometryOriginal(geom.copy()), geometryDeformed(geom),
                                           constraints(std::array<std::vector<double>, 3>{std::vector<double>(),
                                                                                          std::vector<double>(),
                                                                                          std::vector<double>()}),
                                           nConstraints(0), maxIter(30), n(geom.mesh.nVertices()) {
    geom.requireCotanLaplacian();
    lapBeltrami = geom.cotanLaplacian;
    geom.unrequireCotanLaplacian();

    geometryOriginal->requireEdgeCotanWeights();
    geometryOriginal->requireVertexNormals();

    solver = new Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >;

}

ARAP::~ARAP(){
    geometryOriginal->unrequireEdgeCotanWeights();
    geometryOriginal->unrequireVertexNormals();
    delete solver;
}

void ARAP::addConstraint(int id, Vector3 value){
    pairedConstraints.push_back(std::make_pair(id, value));
    nConstraints++;
}

void ARAP::addNormalConstraint(int id, Vector3 value){
    idNormalConstraints.push_back(id);
    normalConstraints.push_back(value);
}

void ARAP::resolveConstraints(){
    sort(pairedConstraints.begin(), pairedConstraints.end(), [](auto& left, auto& right) { return left.first < right.first; });

    idConstraints.resize(pairedConstraints.size());
    constraints[0].resize(pairedConstraints.size());
    constraints[1].resize(pairedConstraints.size());
    constraints[2].resize(pairedConstraints.size());
    for(size_t i = 0; i < pairedConstraints.size(); ++i){
        idConstraints[i] = pairedConstraints[i].first;
        constraints[0][i] = pairedConstraints[i].second.x;
        constraints[1][i] = pairedConstraints[i].second.y;
        constraints[2][i] = pairedConstraints[i].second.z;
    }
}

void ARAP::solve(){

    resolveConstraints();

    std::chrono::time_point<std::chrono::system_clock> startLoop = std::chrono::system_clock::now();
    for(int iter = 0; iter < maxIter; ++iter){
        std::cout << iter << "/" << maxIter << std::endl;

        std::vector<Eigen::Matrix3d> rotations = computeRotations();

        computeNewPositions(rotations);

    }
    std::chrono::time_point<std::chrono::system_clock> endLoop = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endLoop - startLoop;
    std::cout << "compute ARAP : " << elapsed_seconds.count() << std::endl;
}


std::vector<Eigen::Matrix3d> ARAP::computeRotations(){

    auto toEigen = [](Vector3& v) { return Eigen::Vector3d(v.x, v.y, v.z); };

    std::vector<Eigen::Matrix3d> rotations(n);

    for(size_t i = 0; i < n; i++){

        if(std::vector<size_t>::iterator it = std::find(idNormalConstraints.begin(), idNormalConstraints.end(), i); it != idNormalConstraints.end()){
            int id = std::distance(idNormalConstraints.begin(), it);
            Eigen::Vector3d ni = toEigen(geometryOriginal->vertexNormals[i]);
            Eigen::Vector3d ti = toEigen(normalConstraints[id]);
            rotations[i] = findRotationToFitVectors(ni, ti);
            continue;
        }

        Vertex vi = geometryOriginal->mesh.vertex(i);
        Eigen::Matrix3d Si = Eigen::Matrix3d::Zero();

        for(Edge e : vi.adjacentEdges()){
            Vertex vj = e.otherVertex(vi);
            double wij = geometryOriginal->edgeCotanWeights[e];
            Vector3 eij_prime = geometryDeformed.vertexPositions[vi] - geometryDeformed.vertexPositions[vj];
            Vector3 eij = geometryOriginal->vertexPositions[vi] - geometryOriginal->vertexPositions[vj];
            Si += wij * toEigen(eij) * toEigen(eij_prime).transpose();
        }

        // SVD
        Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> svd(Si, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix3d Ri = svd.matrixV() * svd.matrixU().transpose();

        if(Ri.determinant() < 0){
            Eigen::Matrix3d U = svd.matrixU();
            U.coeffRef(0, 2) = -U.coeff(0, 2);
            U.coeffRef(1, 2) = -U.coeff(1, 2);
            U.coeffRef(2, 2) = -U.coeff(2, 2);
            Ri = svd.matrixV() * U.transpose();
        }
        rotations[i] = Ri;
    }
    return rotations;
}

void ARAP::computeNewPositions(std::vector<Eigen::Matrix3d>& R){

    auto toEigen = [](Vector3& v) { return Eigen::Vector3d(v.x, v.y, v.z); };

    Eigen::VectorXd rhsX = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rhsY = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rhsZ = Eigen::VectorXd::Zero(n);

    for(size_t i = 0; i < n; i++){
        Vertex vi = geometryOriginal->mesh.vertex(i);

        if(std::vector<size_t>::iterator it = std::find(idConstraints.begin(), idConstraints.end(), i); it != idConstraints.end()){
            int id = std::distance(idConstraints.begin(), it);
            rhsX[i] = constraints[0][id]; rhsY[i] = constraints[1][id]; rhsZ[i] = constraints[2][id];
            continue;
        }

        for(Edge e : vi.adjacentEdges()){
            Vertex vj = e.otherVertex(vi);
            size_t j = vj.getIndex();
            double wij = geometryOriginal->edgeCotanWeights[e];
            Vector3 pi = geometryOriginal->vertexPositions[vi];
            Vector3 pj = geometryOriginal->vertexPositions[vj];
            Vector3 eij = pi - pj;
            Eigen::Vector3d rotEdge = (R[i] + R[j]) * toEigen(eij);
            rhsX[i] += wij/2 * rotEdge[0];
            rhsY[i] += wij/2 * rotEdge[1];
            rhsZ[i] += wij/2 * rotEdge[2];
        }
    }


    // Solve
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(lapBeltrami, buildConstraintsVector(), true);


    auto solveL = [&decomp, this](Eigen::VectorXd& rhs, int coord) {

        Vector<double> rhsVals1, rhsVals2;
        decomposeVector(decomp, rhs, rhsVals1, rhsVals2);

        Vector<double> combinedRHS = rhsVals1 - decomp.AB * getEigenConstraint(coord);
        Vector<double> Aresult = geometrycentral::solve(decomp.AA, combinedRHS);

        Vector<double> result = reassembleVector(decomp, Aresult, getEigenConstraint(coord));

        return result;

    };

    Vector<double> pX = solveL(rhsX, 0);
    Vector<double> pY = solveL(rhsY, 1);
    Vector<double> pZ = solveL(rhsZ, 2);

    for(size_t i = 0; i < n; ++i){
        geometryDeformed.vertexPositions[i] = Vector3{pX[i], pY[i], pZ[i]};
    }

}

Vector<double> ARAP::getEigenConstraint(int coord){

    Vector<double> eigenConstraints(nConstraints);
    for(int i = 0; i < nConstraints; ++i){
        eigenConstraints[i] = constraints[coord][i];
    }

    return eigenConstraints;
}

Vector<bool> ARAP::buildConstraintsVector(){
    Vector<bool> constraintsVector(n);
    for(size_t i = 0; i < n; ++i)
        constraintsVector[i] = (std::count(idConstraints.begin(), idConstraints.end(), i) == 0);
    return constraintsVector;
}

Eigen::Matrix3d ARAP::findRotationToFitVectors(Eigen::Vector3d a, Eigen::Vector3d b){
    a = a.normalized();
    b = b.normalized();

    Eigen::Vector3d v = a.cross(b);
    double c = a.dot(b);
    Eigen::Matrix3d Vcross;
    Vcross << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + Vcross + Vcross*Vcross*(1.0/(1.0 + c));
    return R;
}
