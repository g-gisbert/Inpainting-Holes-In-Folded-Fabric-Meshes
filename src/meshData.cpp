#include "meshData.hpp"

void MeshInfo::fillInfo() {

}

void MeshInfo::get2Ring() {

    std::vector<size_t> ids(mesh.nVertices());
    for (const Vertex& v: mesh.vertices()) {
        ids.push_back(v.getIndex());
    }

    get2Ring(ids);
}

void MeshInfo::get2Ring(std::vector<size_t>& ids) {

    for (size_t i : ids){
        Vertex v = mesh.vertex(i);
        std::vector<size_t> _1Ring;
        for (const Vertex& v1R : v.adjacentVertices()){
            _1Ring.push_back(v1R.getIndex());
        }
        std::vector<Vertex> current;
        for (const Vertex& v1R : v.adjacentVertices()){
            for (const Vertex& v2R : v1R.adjacentVertices()){
                if (v2R.getIndex() != v.getIndex() &&
                   std::count(_1Ring.begin(),_1Ring.end(),v2R.getIndex()) == 0 &&
                   std::count(current.begin(),current.end(),v2R) == 0){
                    current.push_back(v2R);
                }
            }
        }
        _2Ring[v] = current;
    }
}

void MeshInfo::get2RingEdgeLength() {
    for (const Vertex& v : mesh.vertices()){
        std::vector<double> tmp;
        for (const Vertex& id : _2Ring[v]){
            double value = norm(geom.vertexPositions[v] -
                                geom.vertexPositions[id]);
            tmp.push_back(value);
        }
        _2RingLengths[v] = tmp;
    }
}

void MeshInfo::get1RingEdgeLength() {
    for (const Edge& e : mesh.edges()){
        edgeLengths[e] = norm(geom.vertexPositions[e.firstVertex()] -
                                geom.vertexPositions[e.secondVertex()]);
    }
}

void MeshInfo::getBorder() {
    for (const Vertex& v1 : BL.adjacentVertices()) {
        for (const Vertex& v2 : v1.adjacentVertices()) {
            border.insert(v2);
        }
    }
}

void MeshInfo::getSwallowBorder() {
    for (const Vertex& v : BL.adjacentVertices()) {
        swallowBorder.insert(v);
    }
}

void MeshInfo::getInVertices() {
    for (const Vertex& v : mesh.vertices()) {
        if (border.count(v) == 0) {
            inVertices.insert(v);
        }
    }
}

void MeshInfo::getInEdges() {
    for (const Edge& e : mesh.edges()) {
        if (border.count(e.firstVertex()) == 0 or border.count(e.secondVertex()) == 0) {
            inEdges.insert(e);
        }
    }
}

void MeshInfo::getInFaces() {
    for (const Face& f : mesh.faces()) {
        int n = 0;
        if (border.count(f.halfedge().vertex()))
            n++;
        if (border.count(f.halfedge().next().vertex()))
            n++;
        if (border.count(f.halfedge().next().next().vertex()))
            n++;
        if (n < 3)
            inFaces.insert(f);
    }
}

void MeshInfo::precomputeQ() {
    std::vector<Eigen::Matrix4d> _pQ(mesh.nEdges());

    for (size_t i = 0; i < geom.mesh.nEdges(); ++i) {
        Edge e = mesh.edge(i);
        Face f0 = e.halfedge().face();
        Face f1 = e.halfedge().face();
        double c01 = geom.halfedgeCotanWeight(e.halfedge().next());
        double c02 = geom.halfedgeCotanWeight(e.halfedge().twin().next().next());
        double c03 = geom.halfedgeCotanWeight(e.halfedge().next().next());
        double c04 = geom.halfedgeCotanWeight(e.halfedge().twin().next());

        Eigen::Vector4d K = Eigen::Vector4d {c03+c04, c01+c02, -c01-c03, -c02-c04};
        Eigen::Matrix4d Q = 3.0/(geom.faceArea(f0) + geom.faceArea(f1)) * K * K.transpose();

        _pQ[i] = Q;
    }
    pQ = std::move(_pQ);
}

void MeshInfo::getRepulsiveRadius() {
    for (const Vertex& v : mesh.vertices()) {
        double minEdge = std::numeric_limits<double>::max();
        for (const Edge& e : v.adjacentEdges()) {
            double edgeLength = norm(geom.vertexPositions[e.firstVertex()] - geom.vertexPositions[e.secondVertex()]);
            if (edgeLength < minEdge) {
                minEdge = edgeLength;
            }
        }
        repulsiveRadius[v] = minEdge/1.3;
    }

}


void MeshInfo::precomputeData() {
    getBorder();
    getSwallowBorder();
    getInVertices();
    getInEdges();
    getInFaces();

    get2Ring();
    get2RingEdgeLength();
    get1RingEdgeLength();

    precomputeQ();
    getRepulsiveRadius();
}