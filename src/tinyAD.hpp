#ifndef INC_2DCONTOUR_TINYAD_HPP
#define INC_2DCONTOUR_TINYAD_HPP

#include "pch.h"
#include "meshData.hpp"
#include "utils.hpp"

class TinyADopt {
public:


    TinyADopt(MeshInfo& md) : m_alpha(0.1), m_beta(0.0001), max_iters(1000), convergence_eps(1e-2), ref_md(md) {

    }

    void run();

    double m_alpha;
    double m_beta;

private:
    int max_iters;
    double convergence_eps;
    MeshInfo& ref_md;

    TinyAD::LinearSolver<double, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> solver;
};

#endif //INC_2DCONTOUR_TINYAD_HPP
