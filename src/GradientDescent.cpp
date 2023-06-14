#include "GradientDescent.hpp"



void GradientDescent::step(const bool bendingEnergy, const bool enableCollisions) noexcept {
    std::fill(grad.begin(), grad.end(), Vector3::zero());

    VertexData<Vector3> &vp = meshData.geom.vertexPositions;

    // Edge Length
    double Elength = 0.0;
    for (const Edge &e: meshData.mesh.edges()) {
        size_t id1 = e.firstVertex().getIndex();
        size_t id2 = e.secondVertex().getIndex();
        double l = meshData.edgeLengths[e];
        Vector3 u = -vp[id1] + vp[id2];
        double normU = norm(u);
        u /= norm(u);
        grad[id1] += alpha * (l - normU) * u;
        grad[id2] -= alpha * (l - normU) * u;
        Elength += (l - normU) * (l - normU);
    }
    std::cout << "Elength : " << Elength << std::endl;

    // 2 Ring
    if (!bendingEnergy) {
        double Ebending = 0.0;
        for (const Vertex &v1: meshData.mesh.vertices()) {
            for (size_t i = 0; i < meshData._2Ring[v1].size(); ++i) {
                Vertex v2 = meshData._2Ring[v1][i];
                size_t id1 = v1.getIndex();
                size_t id2 = v2.getIndex();
                double l = meshData._2RingLengths[v1][i];
                Vector3 u = -vp[id1] + vp[id2];
                double normU = norm(u);
                u /= norm(u);
                grad[id1] += beta * (l - normU) * u;
                grad[id2] -= beta * (l - normU) * u;
            }
        }
        std::cout << "Ebending : " << Ebending << std::endl;
    }
    else {
        // Bending Energy
        double Ebending = 0.0;
        for (const Edge &e: meshData.mesh.edges()) {

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
            grad[x0] += 0.001 * beta * Vector3{G(0, 0), G(0, 1), G(0, 2)};
            grad[x1] += 0.001 * beta * Vector3{G(1, 0), G(1, 1), G(1, 2)};
            grad[x2] += 0.001 * beta * Vector3{G(2, 0), G(2, 1), G(2, 2)};
            grad[x3] += 0.001 * beta * Vector3{G(3, 0), G(3, 1), G(3, 2)};
            Ebending += (0.5 * x.transpose() * meshData.pQ[e.getIndex()] * x).trace();
        }
        std::cout << "Ebending : " << Ebending << std::endl;
    }

    if (enableCollisions) {
        updateSpatialHash();
        detectCollisions();
    }

    update();
}

void GradientDescent::update() noexcept {
    VertexData<Vector3>& vp = meshData.geom.vertexPositions;
    for (const Vertex& v : meshData.inVertices) {
        vp[v] = vp[v] - timeStep*grad[v.getIndex()];
    }
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

    // Self-intersections
    for (const auto& [hashValue, vertexList] : spatialHashing) {
        for (const Point3d& p1 : vertexList) {
            for (const Point3d& p2 : vertexList) {
                if (p1 == p2)
                    continue;
                if (norm(p1.v - p2.v) < meshData.repulsiveRadius[p1.id]) {
                    Vector3 u = (p1.v - p2.v).normalize();
                    grad[p1.id] -= 0.0001*(mu / norm2(p1.v - p2.v)) * u;
                }

            }
        }
    }

}

void GradientDescent::updateSpatialHash() noexcept {
    VertexData<Vector3>& vp = meshData.geom.vertexPositions;
    spatialHashing.clear();

    for (const Vertex& v : meshData.mesh.vertices()) {
        int h = hash(vp[v]);
        spatialHashing[h].push_back(Point3d{vp[v], v.getIndex()});
    }

}

int GradientDescent::hash(const Vector3 &v) const noexcept {
    int x = int(v.x / gridSize);
    int y = int(v.y / gridSize);
    int z = int(v.z / gridSize);
    return x + y * 10 + z * 100;
}

void GradientDescent::parseObjects() {
    std::ifstream file("../config.txt");
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "sphere") {
            double x, y, z, r;
            iss >> x >> y >> z >> r;
            objects.push_back(std::make_unique<Sphere>(Vector3{x, y, z}, r));
        } else if (type == "cylinder") {
            double x, y, z, r, h, a1, a2, a3;
            iss >> x >> y >> z >> r >> h >> a1 >> a2 >> a3;
            objects.push_back(std::make_unique<Cylinder>(Vector3{x, y, z}, r, h, Utils::eulerAngles(a1, a2, a3)));
        } else if (type[0] != '#') {
            std::cout << "Invalid type detected: " << type << std::endl;
        }
    }
}