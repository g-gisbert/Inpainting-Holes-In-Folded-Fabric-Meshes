#include "optimizer.hpp"

int graddescent(Optimizer::Func f, Eigen::VectorXd& x) {
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(2 * f.chain.size());
    int i = 1;
    for(int j = 0; j < 100000; ++j) {
        f(x, grad);
        if (j%1000 == 0)
            Chain::toSVG(x, grad, i++);
        x = x - 0.1*grad;

    }

    return 30000;
}

void Optimizer::lbfgs() {
    Func func(ref_chain);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(2 * ref_chain.size());

    for (int i = 0; i < int(ref_chain.size()); ++i) {
        x[2*i+0] = ref_chain.getPosition(i).x;
        x[2*i+1] = ref_chain.getPosition(i).y;
    }
    //std::cout << x;
    double E = 0.0;

    auto start = std::chrono::high_resolution_clock::now();
    //int nIter = graddescent(func, x);
    int nIter = m_solver->minimize(func, x, E);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "nIter = " << nIter << ", E = " << E << std::endl;
    auto elapsed_time = duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Elapsed time: " << elapsed_time << " microseconds" << std::endl;

    for (int i = 0; i < int(ref_chain.size()); ++i) {
        ref_chain.getPosition(i) = Vector2{x[2*i], x[2*i+1]};
    }
}

double Optimizer::Func::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    grad = Eigen::VectorXd::Zero(2 * chain.size());

    constexpr double alpha = 0.001;
    constexpr double beta = 0.01;

    const int N = int(chain.size());

    // Edge Lengths gradient
   double Eiso = 0.0f;
    for (int i = 0; i < N; ++i) {
        double l = chain.getEdgeLength(i);
        Vector2 posV1 = Vector2{x[2*i], x[2*i+1]};
        Vector2 posV2 = Vector2{x[2*((i+1)%N)], x[2*((i+1)%N)+1]};
        Vector2 u = (posV2 - posV1).normalize();

        Vector2 vec = 2 * alpha * ( l - norm(posV1 - posV2)) * u;

        grad[2 * i + 0] += vec.x;
        grad[2 * i + 1] += vec.y;

        grad[2 * ((i+1)%N) + 0] -= vec.x;
        grad[2 * ((i+1)%N) + 1] -= vec.y;

        double edgeLength = norm(posV1 - posV2);
        Eiso += (edgeLength-l)*(edgeLength-l);
    }

    //std::cout << "Eiso : " << alpha*Eiso << std::endl;

    // Angles Gradient
    double Eangle = 0.0f;
    for (int i = 0; i < N; ++i){
        int id = ((i-1) < 0) ? i-1+N : i-1;
        int id2 = (i+1)%N;

        Vector2 posV1 = Vector2{x[2*id], x[2*id+1]};
        Vector2 posV2 = Vector2{x[2*i], x[2*i+1]};
        Vector2 posV3 = Vector2{x[2*id2], x[2*id2+1]};

        Vector2 incomingVector = (posV2 - posV1).normalize();
        Vector2 outgoingVector = (posV3 - posV2).normalize();

        double angle = atan2(-incomingVector.y, -incomingVector.x) - atan2(outgoingVector.y, outgoingVector.x);
        if (angle < 0) { angle += 2 * M_PI; }


        double debugl = chain.getAngle(i);
        double front = beta * 2 * (chain.getAngle(i) - angle);
        Vector2 cross1 = Vector2{incomingVector.y, -incomingVector.x} / norm2(incomingVector);
        Vector2 cross2 = Vector2{outgoingVector.y, -outgoingVector.x} / norm2(outgoingVector);

        grad[2*i+0] -= front * (-cross1.x - cross2.x);
        grad[2*i+1] -= front * (-cross1.y - cross2.y);

        grad[2 * id + 0] -= front * cross1.x;
        grad[2 * id + 1] -= front * cross1.y;

        grad[2 * id2 + 0] -= front * cross2.x;
        grad[2 * id2 + 1] -= front * cross2.y;

        Eangle += (chain.getAngle(i) - angle)*(chain.getAngle(i) - angle);
    }

    grad[0] = 0;
    grad[1] = 0;

    //std::cout << "max c : " << maxC << " || max d : " << maxD << " || Eiso : " << opt.m_alpha*Eiso << " || Eangle : " << opt.m_beta*Eangle << std::endl;
    std::cout << alpha*Eiso << " || Eangle : " << beta*Eangle << std::endl;
    double e = alpha*Eiso + beta*Eangle;
    return e;
}