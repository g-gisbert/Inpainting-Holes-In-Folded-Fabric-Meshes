#include "utils.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


namespace Utils {

    double averageBlEdgeLength(const BoundaryLoop& BL, const VertexPositionGeometry& geom) {
        int nEdges = 0;
        double average = 0.0;
        for (const Edge& e : BL.adjacentEdges()) {
            average += norm(geom.vertexPositions[e.firstVertex()] - geom.vertexPositions[e.secondVertex()]);
            ++nEdges;
        }
        average /= nEdges;
        return average;
    }


    void ARAP(VertexPositionGeometry& geom, VertexPositionGeometry& geomValues,
              const MeshInfo& meshData, const std::unordered_map<size_t, size_t>& matching, const double c) {

        geomValues.requireVertexNormals();

        for (const Vertex& v : geom.mesh.vertices()) {
            geom.vertexPositions[v] = geom.vertexPositions[v] * c;
        }

        class ARAP arapHolder(geom);

        for (const Vertex& v : meshData.border) {
            size_t id2D = v.getIndex();
            size_t id3D = matching.at(id2D);
            arapHolder.addConstraint(id2D, geomValues.vertexPositions[id3D]);
        }

        geomValues.unrequireVertexNormals();
        arapHolder.solve();
    }

    AABB getAABB(const VertexPositionGeometry& geom) {

        const VertexData<Vector3>& vp = geom.vertexPositions;
        Vector3 min = vp[0];
        Vector3 max = vp[0];

        for (const Vertex& v : geom.mesh.vertices()) {
            const Vector3 pos = vp[v];
            if (pos.x < min.x)
                min.x = pos.x;
            if (pos.y < min.y)
                min.y = pos.y;
            if (pos.z < min.z)
                min.z = pos.z;
            if (pos.x > max.x)
                max.x = pos.x;
            if (pos.y > max.y)
                max.y = pos.y;
            if (pos.z > max.z)
                max.z = pos.z;
        }
        return AABB{min, max};
    }

    std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
        createCylinder(double radius, double height, int N, Vector3 center, Eigen::Matrix3d R) {

        std::vector<Vector3> pts(2*N);
        for (int i = 0; i < N; ++i) {
            double angle = 2.0 * i * M_PI / N;
            double x = radius * cos(angle);
            double y = radius * sin(angle);
            pts[2*i]   = Vector3{x, height/2, y};
            pts[2*i+1] = Vector3{x, -height/2, y};
        }

        std::vector<std::vector<size_t>> faces;
        for (size_t i = 0; i < size_t(N); ++i) {
            size_t next = (i + 1) % N;

            // Side faces
            std::vector<size_t> face = {i*2, next*2, i*2+1};
            faces.push_back(face);

            face = {next*2, next*2+1, i*2+1};
            faces.push_back(face);

            if (i >= size_t(N) - 2)
                continue;
            face = {0, 2*(next+1), 2*next};
            faces.push_back(face);
            face = {1, 2*next+1, 2*(next+1)+1};
            faces.push_back(face);
        }

        for (int i = 0; i < 2*N; ++i) {
            Eigen::Vector3d tmp = R * Eigen::Vector3d{pts[i].x, pts[i].y, pts[i].z};
            pts[i] = Vector3{tmp[0], tmp[1], tmp[2]} + center;
        }

        return makeManifoldSurfaceMeshAndGeometry(faces, pts);
    }

    std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
        createIcoSphere(double scale, int level, Vector3 center) {

        const float X=.525731112119133606f;
        const float Z=.850650808352039932f;
        const float N=0.f;

        std::vector<Vector3> pts = {Vector3{-X,N,Z}, Vector3{X,N,Z}, Vector3{-X,N,-Z}, Vector3{X,N,-Z},
                                    Vector3{N,Z,X}, Vector3{N,Z,-X}, Vector3{N,-Z,X}, Vector3{N,-Z,-X},
                                    Vector3{Z,X,N}, Vector3{-Z,X, N}, Vector3{Z,-X,N}, Vector3{-Z,-X, N}};

        std::vector<std::vector<size_t>> faces = {{0,1,4},{0,4,9},{9,4,5},{4,8,5},{4,1,8},
                                                  {8,1,10},{8,10,3},{5,8,3},{5,3,2},{2,3,7},
                                                  {7,3,10},{7,10,6},{7,6,11},{11,6,0},{0,6,1},
                                                  {6,10,1},{9,11,0},{9,2,11},{9,5,2},{7,11,2}};

        std::unique_ptr<ManifoldSurfaceMesh> meshIco;
        std::unique_ptr<VertexPositionGeometry> geometryIco;
        std::tie(meshIco, geometryIco) = makeManifoldSurfaceMeshAndGeometry(faces, pts);

        for(int lvl = 0; lvl < level; ++lvl){
            for(Edge e : meshIco->edges()){
                Vertex v0 = e.firstVertex();
                Vertex v1 = e.secondVertex();
                Vector3 vPos0 = geometryIco->inputVertexPositions[v0];
                Vector3 vPos1 = geometryIco->inputVertexPositions[v1];

                Halfedge he = meshIco->insertVertexAlongEdge(e);

                geometryIco->inputVertexPositions[he.vertex()] = (vPos0 + vPos1) * 0.5f / norm((vPos0 + vPos1) * 0.5f);
            }

            for(Face f : meshIco->faces()){
                Halfedge he = f.halfedge();
                if(he.vertex().degree() != 2){
                    he = he.next();
                }
                meshIco->connectVertices(he, he.next().next());
                meshIco->connectVertices(he.next().next().twin().next(), he.next().next().twin().next().next().next());
                meshIco->connectVertices(he.next().next().twin().next().next(), he.next().next().twin().next().next().next().next());
            }
        }

        for (const Vertex& v : meshIco->vertices()) {
            geometryIco->vertexPositions[v] = geometryIco->vertexPositions[v]*scale + center;
        }

        return {std::move(meshIco), std::move(geometryIco)};
    }


    Eigen::Matrix3d eulerAngles(double psi, double phi, double theta) {
        return Eigen::Matrix3d{{cos(psi)*cos(phi)-sin(psi)*cos(theta)*sin(phi), -cos(psi)*sin(phi)-sin(psi)*cos(theta)*cos(phi), sin(psi)*sin(theta)},
                               {sin(psi)*cos(phi)+cos(psi)*cos(theta)*sin(phi), -sin(psi)*sin(phi)+cos(psi)*cos(theta)*cos(phi), -cos(psi)*sin(theta)},
                               {sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta)}};
    }


    Eigen::Vector3d to_eigen(const Vector3& v) {
        return Eigen::Vector3d{v.x, v.y, v.z};
    }

    Vector3 to_geometrycentral(const Eigen::Vector3d& v) {
        return Vector3{v[0], v[1], v[2]};
    }

}