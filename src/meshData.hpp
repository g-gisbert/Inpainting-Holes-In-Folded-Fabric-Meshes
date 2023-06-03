#ifndef INC_2DCONTOUR_MESHDATA_HPP
#define INC_2DCONTOUR_MESHDATA_HPP

#include "pch.h"

class MeshInfo {
public:
    MeshInfo(ManifoldSurfaceMesh& _mesh, VertexPositionGeometry& _geom) : mesh(_mesh), geom(_geom),
        BL(mesh.boundaryLoop(0)), _2Ring(_mesh), edgeLengths(_mesh), _2RingLengths(_mesh), repulsiveRadius(_mesh) {};

    void fillInfo();

    // Mesh info
    ManifoldSurfaceMesh& mesh;
    VertexPositionGeometry& geom;

    // selection data
    BoundaryLoop BL;
    std::set<Vertex> swallowBorder;
    std::set<Vertex> border; // ok (thick)
    std::set<Vertex> inVertices;
    std::set<Edge> inEdges;
    std::set<Face> inFaces; // ok
    VertexData<std::vector<Vertex>> _2Ring; // ok

    // quantitative data
    EdgeData<double> edgeLengths; // ok
    VertexData<std::vector<double>> _2RingLengths; // ok
    std::vector<Eigen::Matrix4d> pQ; // ok
    VertexData<double> repulsiveRadius;

    void precomputeData();

private:

    void getBorder();
    void getSwallowBorder();
    void getInVertices();
    void getInEdges();
    void getInFaces();

    void get2Ring();
    void get2Ring(std::vector<size_t>& ids);
    void get2RingEdgeLength();
    void get1RingEdgeLength();

    void precomputeQ();

    void getRepulsiveRadius();

};


#endif //INC_2DCONTOUR_MESHDATA_HPP
