#include "GradientDescent.hpp"



void GradientDescent::step() noexcept {
    std::fill(grad.begin(), grad.end(), Vector3::zero());

    VertexData<Vector3>& vp = meshData.geom.vertexPositions;

    //std::ofstream file("data.txt", std::ios::app);


    // Edge Length
    double Elength = 0.0;
    for (const Edge& e : meshData.mesh.edges()) {
        size_t id1 = e.firstVertex().getIndex();
        size_t id2 = e.secondVertex().getIndex();
        double l = lambda*meshData.edgeLengths[e];
        Vector3 u = -vp[id1] + vp[id2];
        double normU = norm(u);
        u /= norm(u);
        grad[id1] += alpha*(l - normU) * u;
        grad[id2] -= alpha*(l - normU) * u;
        Elength += (l - normU)*(l - normU);
    }
    std::cout << "Elength : " << Elength << std::endl;

    // 2 Ring
    for (const Vertex& v1 : meshData.mesh.vertices()) {
        for (size_t i = 0; i < meshData._2Ring[v1].size(); ++i) {
            Vertex v2 = meshData._2Ring[v1][i];
            size_t id1 = v1.getIndex();
            size_t id2 = v2.getIndex();
            double l = meshData._2RingLengths[v1][i];
            Vector3 u = -vp[id1] + vp[id2];
            double normU = norm(u);
            u /= norm(u);
            grad[id1] += beta*(l - normU) * u;
            grad[id2] -= beta*(l - normU) * u;
        }
    }

    // Bending Energy
    double Ebending = 0.0;
    for (const Edge& e : meshData.mesh.edges()) {

        size_t x0 = e.halfedge().vertex().getIndex();
        size_t x1 = e.halfedge().next().vertex().getIndex();
        size_t x2 = e.halfedge().next().next().vertex().getIndex();
        size_t x3 = e.halfedge().twin().next().next().vertex().getIndex();
        Vector3 x0p = vp[x0];
        Vector3 x1p = vp[x1];
        Vector3 x2p = vp[x2];
        Vector3 x3p = vp[x3];
        Eigen::Matrix<double, 4, 3> x;
        x << x0p.x, x0p.y, x0p.z, x1p.x, x1p.y, x1p.z, x2p.x, x2p.y, x2p.z, x3p.x, x3p.y, x3p.z;

        Eigen::Matrix<double, 4, 3> G = 2 * meshData.pQ[e.getIndex()] * x;
        grad[x0] += 0.001*gamma * Vector3{G(0, 0), G(0, 1), G(0, 2)};
        grad[x1] += 0.001*gamma * Vector3{G(1, 0), G(1, 1), G(1, 2)};
        grad[x2] += 0.001*gamma * Vector3{G(2, 0), G(2, 1), G(2, 2)};
        grad[x3] += 0.001*gamma * Vector3{G(3, 0), G(3, 1), G(3, 2)};
        Ebending += (0.5 * x.transpose() * meshData.pQ[e.getIndex()] * x).trace();
    }
    std::cout << "Ebending : " << Ebending << std::endl;

    //updateSpatialHash();
    //detectCollisions();
    /*std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "detectCollisions() : " << elapsed_seconds.count() << std::endl;*/
    //detectCollisions();
    update();

    /*meshData.geom.requireVertexGaussianCurvatures();
    meshData.geom.requireVertexMeanCurvatures();
    double average = 0;
    double worst = 0.0;
    double variance = 0.0;
    for (const Vertex& v : meshData.mesh.vertices()) {
        average += abs(meshData.geom.vertexMeanCurvatures[v]);
        if (abs(meshData.geom.vertexMeanCurvatures[v]) > worst)
            worst = abs(meshData.geom.vertexMeanCurvatures[v]);
    }
    average /= meshData.mesh.nVertices();

    for (const Vertex& v : meshData.mesh.vertices()) {
        variance += (abs(meshData.geom.vertexMeanCurvatures[v]) - average)*(abs(meshData.geom.vertexMeanCurvatures[v]) - average);
    }
    variance /= meshData.mesh.nVertices();
    file << std::to_string(iter++) << " " << average << " " << variance << std::endl;
    meshData.geom.unrequireVertexGaussianCurvatures();
    file.close();*/
}

void GradientDescent::update() noexcept {
    VertexData<Vector3>& vp = meshData.geom.vertexPositions;
    for (const Vertex& v : meshData.inVertices) {
        vp[v] = vp[v] - timeStep*grad[v.getIndex()];
    }
    //polyscope::getSurfaceMesh("ARAP")->updateVertexPositions(meshData.geom.vertexPositions);
}

void GradientDescent::detectCollisions() noexcept {

    VertexData<Vector3>& vp = meshData.geom.vertexPositions;

    // Cloth object collision
    for (std::unique_ptr<Shape>& sp : objects) {
        for (const Vertex& v : meshData.mesh.vertices()) {
            if (sp->isInside(vp[v])){
                vp[v] = sp->closestSurfacePoint(vp[v]);
            }
        }
    }

    // TODO : Vector3 & id struct + movable sphere
    for (const auto& [hashValue, vertexList] : spatialHashing) {
        for (const Point3d& p1 : vertexList) {
            for (const Point3d& p2 : vertexList) {
                if (p1 == p2)
                    continue;
                if (norm(p1.v - p2.v) < meshData.repulsiveRadius[p1.id]) {
                    std::cout << "intersection : " << p1.id << ", " << p2.id << " / " << norm(p1.v - p2.v) << std::endl;
                    Vector3 u = (p1.v - p2.v).normalize();
                    grad[p1.id] -= 0.001*(mu / norm2(p1.v - p2.v)) * u;
                }

            }
        }
    }


}

void GradientDescent::updateSpatialHash() noexcept {
    VertexData<Vector3>& vp = meshData.geom.vertexPositions;
    spatialHashing.clear();
    /*for (const Face& f : meshData.mesh.faces()) {
        Halfedge he = f.halfedge();
        Vertex v1 = he.vertex();
        Vertex v2 = he.next().vertex();
        Vertex v3 = he.next().next().vertex();
        int hash1 = hash(vp[v1]);
        int hash2 = hash(vp[v2]);
        int hash3 = hash(vp[v3]);
        Triangle t{vp[v1], vp[v2], vp[v3], f.getIndex()};
        spatialHashing[hash1].insert(t);
        spatialHashing[hash2].insert(t);
        spatialHashing[hash3].insert(t);
    }*/

    for (const Vertex& v : meshData.mesh.vertices()) {
        int h = hash(vp[v]);
        spatialHashing[h].push_back(Point3d{vp[v], v.getIndex()});
    }

    /*VertexData<int> spatial(meshData.mesh);
    for (const Vertex& v : meshData.mesh.vertices()) {
        spatial[v] = hash(vp[v]);
    }
    polyscope::getSurfaceMesh("ARAP")->addVertexScalarQuantity("hash", spatial);*/

}

int GradientDescent::hash(const Vector3 &v) const noexcept {
    int x = int(v.x / gridSize);
    int y = int(v.y / gridSize);
    int z = int(v.z / gridSize);
    return x + y * 10 + z * 100;
}
