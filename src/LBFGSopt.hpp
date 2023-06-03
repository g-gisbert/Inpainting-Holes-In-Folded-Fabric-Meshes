#ifndef INC_2DCONTOUR_LBFGSOPT_HPP
#define INC_2DCONTOUR_LBFGSOPT_HPP

#include <LBFGS.h>
#include "meshData.hpp"
#include "pch.h"

using namespace LBFGSpp;

class LBFGSopt {
public:

    struct Func {
        MeshInfo& md;
        Func(MeshInfo& md_) : md(md_) {}
        double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
    };

    LBFGSopt(MeshInfo& md) : m_alpha(0.1), m_beta(0.0001), ref_md(md) {
        m_param.epsilon = 1e-8;
        m_param.epsilon_rel = 1e-8;
        m_param.max_linesearch = 20;
        m_param.max_iterations = 500;
        //m_param.linesearch = 3;
        m_param.wolfe = 0.9;
        //m_param.ftol = 1e-6;
        m_solver = std::make_unique<LBFGSSolver<double>>(m_param);
    }

    void run();

    double m_alpha;
    double m_beta;

private:
    MeshInfo& ref_md;

    LBFGSParam<double> m_param;
    std::unique_ptr<LBFGSSolver<double>> m_solver;
};


#endif //INC_2DCONTOUR_LBFGSOPT_HPP
