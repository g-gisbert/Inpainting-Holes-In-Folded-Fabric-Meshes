#ifndef INC_2DCONTOUR_OPTIMIZER_HPP
#define INC_2DCONTOUR_OPTIMIZER_HPP

#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>
#include "chain.hpp"
#include <chrono>

using namespace LBFGSpp;

class Optimizer {
public:

    class Func {
    public:
        Chain& chain;
    public:
        Func(Chain& chain_) : chain(chain_) {}
        double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
    };


    Optimizer(Chain& chain) : ref_chain(chain), m_alpha(0.1), m_beta(0.1) { //alpha=0.1, beta=0.001
        m_param.epsilon = 1e-8;
        m_param.epsilon_rel = 1e-8;
        m_param.max_linesearch = 20;
        m_param.max_iterations = 1000;
        //m_param.linesearch = 3;
        m_param.wolfe = 0.9;
        //m_param.ftol = 1e-6;
        m_solver = std::make_unique<LBFGSSolver<double>>(m_param);
    }

    void setAlphaBeta(double a, double b) {
        m_alpha = a;
        m_beta = b;
    }

    void lbfgs();

private:
    Chain& ref_chain;

    double m_alpha;
    double m_beta;

    LBFGSParam<double> m_param;
    std::unique_ptr<LBFGSSolver<double>> m_solver;
};

#endif //INC_2DCONTOUR_OPTIMIZER_HPP
