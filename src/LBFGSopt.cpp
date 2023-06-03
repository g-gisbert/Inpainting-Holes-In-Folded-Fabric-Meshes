#include "LBFGSopt.hpp"

void LBFGSopt::run() {
    Func func(ref_md);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(3 * ref_md.mesh.nVertices());

    for (int i = 0; i < int(ref_md.mesh.nVertices()); ++i) {
        x[3*i+0] = ref_md.geom.vertexPositions[i].x;
        x[3*i+1] = ref_md.geom.vertexPositions[i].y;
        x[3*i+2] = ref_md.geom.vertexPositions[i].z;
    }
    //std::cout << x;
    double E = 0.0;

    auto start = std::chrono::high_resolution_clock::now();
    int nIter = m_solver->minimize(func, x, E);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "nIter = " << nIter << ", E = " << E << std::endl;
    auto elapsed_time = duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Elapsed time: " << elapsed_time << " microseconds" << std::endl;

    for (int i = 0; i < int(ref_md.mesh.nVertices()); ++i) {
        ref_md.geom.vertexPositions[i] = Vector3{x[3*i], x[3*i+1], x[3*i+2]};
    }
}

double LBFGSopt::Func::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    grad = Eigen::VectorXd::Zero(3 * md.mesh.nVertices());

    constexpr double alpha = 0.1;
    constexpr double beta = 0.01;//0.01;

    const int N = int(md.mesh.nVertices());

    // Edge Lengths gradient
    double Eiso = 0.0f;

    for (const Edge& e : md.mesh.edges()) {
        size_t id1 = e.firstVertex().getIndex();
        size_t id2 = e.secondVertex().getIndex();
        double l = md.edgeLengths[e];
        Vector3 posV1 = Vector3{x[3*id1+0], x[3*id1+1], x[3*id1+2]};
        Vector3 posV2 = Vector3{x[3*id2+0], x[3*id2+1], x[3*id2+2]};
        Vector3 u = -posV1 + posV2;
        double normU = norm(u);
        u /= norm(u);
        grad[3*id1+0] += (alpha*(l - normU) * u).x;
        grad[3*id1+1] += (alpha*(l - normU) * u).y;
        grad[3*id1+2] += (alpha*(l - normU) * u).z;

        grad[3*id2+0] -= (alpha*(l - normU) * u).x;
        grad[3*id2+1] -= (alpha*(l - normU) * u).y;
        grad[3*id2+2] -= (alpha*(l - normU) * u).z;

        double edgeLength = norm(posV1 - posV2);
        Eiso += (edgeLength-l)*(edgeLength-l);
    }

    //std::cout << "Eiso : " << alpha*Eiso << std::endl;

    // Bending Energy
    double Eangle = 0;
    for (const Edge& e : md.mesh.edges()) {

        size_t x0 = e.halfedge().vertex().getIndex();
        size_t x1 = e.halfedge().next().vertex().getIndex();
        size_t x2 = e.halfedge().next().next().vertex().getIndex();
        size_t x3 = e.halfedge().twin().next().next().vertex().getIndex();
        Vector3 x0p = Vector3{x[3*x0+0], x[3*x0+1], x[3*x0+2]};
        Vector3 x1p = Vector3{x[3*x1+0], x[3*x1+1], x[3*x1+2]};;
        Vector3 x2p = Vector3{x[3*x2+0], x[3*x2+1], x[3*x2+2]};;
        Vector3 x3p = Vector3{x[3*x3+0], x[3*x3+1], x[3*x3+2]};;
        Eigen::Matrix<double, 4, 3> xalt;
        xalt << x0p.x, x0p.y, x0p.z, x1p.x, x1p.y, x1p.z, x2p.x, x2p.y, x2p.z, x3p.x, x3p.y, x3p.z;

        Eigen::Matrix<double, 4, 3> G = 2 * md.pQ[e.getIndex()] * xalt;
        grad[3*x0+0] += (0.001*beta * Vector3{G(0, 0), G(0, 1), G(0, 2)}).x;
        grad[3*x0+1] += (0.001*beta * Vector3{G(0, 0), G(0, 1), G(0, 2)}).y;
        grad[3*x0+2] += (0.001*beta * Vector3{G(0, 0), G(0, 1), G(0, 2)}).z;

        grad[3*x1+0] += (0.001*beta * Vector3{G(1, 0), G(1, 1), G(1, 2)}).x;
        grad[3*x1+1] += (0.001*beta * Vector3{G(1, 0), G(1, 1), G(1, 2)}).y;
        grad[3*x1+2] += (0.001*beta * Vector3{G(1, 0), G(1, 1), G(1, 2)}).z;

        grad[3*x2+0] += (0.001*beta * Vector3{G(2, 0), G(2, 1), G(2, 2)}).x;
        grad[3*x2+1] += (0.001*beta * Vector3{G(2, 0), G(2, 1), G(2, 2)}).y;
        grad[3*x2+2] += (0.001*beta * Vector3{G(2, 0), G(2, 1), G(2, 2)}).z;

        grad[3*x3+0] += (0.001*beta * Vector3{G(3, 0), G(3, 1), G(3, 2)}).x;
        grad[3*x3+1] += (0.001*beta * Vector3{G(3, 0), G(3, 1), G(3, 2)}).y;
        grad[3*x3+2] += (0.001*beta * Vector3{G(3, 0), G(3, 1), G(3, 2)}).z;

        Eangle += (0.5 * xalt.transpose() * md.pQ[e.getIndex()] * xalt).trace();
    }

    // TODO : zero grad thick border
    for (const Vertex& v : md.border) {
        size_t id = v.getIndex();
        grad[3*id+0] = 0;
        grad[3*id+1] = 0;
        grad[3*id+2] = 0;
    }

    //std::cout << "max c : " << maxC << " || max d : " << maxD << " || Eiso : " << opt.m_alpha*Eiso << " || Eangle : " << opt.m_beta*Eangle << std::endl;
    std::cout << "Eiso : " << alpha*Eiso << " || Eangle : " << 0.001*beta*Eangle << std::endl;
    double e = alpha*Eiso + 0.001*beta*Eangle;
    return e;
}