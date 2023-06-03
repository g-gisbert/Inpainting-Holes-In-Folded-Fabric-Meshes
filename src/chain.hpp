#ifndef INC_2DCONTOUR_CHAIN_HPP
#define INC_2DCONTOUR_CHAIN_HPP

#include "pch.h"
#include "plane.hpp"
#include "initialization.hpp"
#include "simple_svg_1.0.0.hpp"

class Chain {
public:

    void readFromMesh(std::string filename);

    void setPositions(Vector2 p);
    void saveTargetEdgeLengths();
    void saveTargetAngles();

    Vector2& getPosition(int i) {
        return positions[i];
    }

    double getEdgeLength(int i) {
        return edgeLengths[i];
    }

    double getAngle(int i) {
        return angles[i];
    }

    void printAngles() const {
        for (size_t i = 0; i < size()-2; ++i) {
            std::cout << angles[i] << std::endl;
        }
    }

    size_t size() const {
        return positions.size();
    }

    void setPositionsFromMesh(const std::string& filename);

    void computeSVGTransform();
    void toSVG();
    static void toSVG(Eigen::VectorXd& x, Eigen::VectorXd& grad, int k);

private:
    std::vector<Vector2> positions;
    std::vector<double> edgeLengths;
    std::vector<double> angles;
    int frame{};

    // Mesh info
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;

    //scale + translation
    static double scale;
    static Vector2 t;

};


#endif //INC_2DCONTOUR_CHAIN_HPP
