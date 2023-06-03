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

    void ARAP(class ARAP& arapHolder, VertexPositionGeometry& geom, VertexPositionGeometry& geomValues,
          const MeshInfo& meshData, const std::unordered_map<size_t, size_t>& matching) {

        for (const Vertex& v : meshData.border) {
            size_t id2D = v.getIndex();
            size_t id3D = matching.at(id2D);
            arapHolder.addConstraint(id2D, geomValues.vertexPositions[id3D]);
        }

        arapHolder.solve();
    }

    void harmonicSurface(VertexPositionGeometry& geom, VertexPositionGeometry& geomValues,
                         const MeshInfo& meshData, const std::unordered_map<size_t, size_t>& matching) {

        geom.requireCotanLaplacian();
        Eigen:SparseMatrix<double> L = geom.cotanLaplacian;

        Vector<double> bcValsX = Vector<double>::Zero(meshData.border.size());
        Vector<double> bcValsY = Vector<double>::Zero(meshData.border.size());
        Vector<double> bcValsZ = Vector<double>::Zero(meshData.border.size());
        Vector<bool> isInterior = Vector<bool>::Constant(geom.mesh.nVertices(), true);

        int j = 0;
        for (const Vertex& v : geom.mesh.vertices()) {
            if (meshData.border.count(v) > 0) {
                isInterior[v.getIndex()] = false;
                size_t id2D = v.getIndex();
                size_t id3D = matching.at(id2D);
                bcValsX[j] = geomValues.vertexPositions[id3D].x;
                bcValsY[j] = geomValues.vertexPositions[id3D].y;
                bcValsZ[j] = geomValues.vertexPositions[id3D].z;
                j++;
            }
        }
        BlockDecompositionResult<double> decomp = blockDecomposeSquare(L, isInterior, true);

        Vector<double> combinedRHSX = -decomp.AB * bcValsX;
        Vector<double> combinedRHSY = -decomp.AB * bcValsY;
        Vector<double> combinedRHSZ = -decomp.AB * bcValsZ;

        Vector<double> AresultX = solve(decomp.AA, combinedRHSX);
        Vector<double> resultX = reassembleVector(decomp, AresultX, bcValsX);
        Vector<double> AresultY = solve(decomp.AA, combinedRHSY);
        Vector<double> resultY = reassembleVector(decomp, AresultY, bcValsY);
        Vector<double> AresultZ = solve(decomp.AA, combinedRHSZ);
        Vector<double> resultZ = reassembleVector(decomp, AresultZ, bcValsZ);

        for (const Vertex& v : geom.mesh.vertices()) {
            geom.vertexPositions[v] = Vector3{resultX[v.getIndex()], resultY[v.getIndex()], resultZ[v.getIndex()]};
        }

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

    bool checkIntersection(const Triangle& t1, const Triangle& t2) noexcept {
        return CGAL::do_intersect(K::Triangle_3{K::Point_3{t1.v1.x, t1.v1.y, t1.v1.z},K::Point_3{t1.v2.x, t1.v2.y, t1.v2.z},
                                         K::Point_3{t1.v3.x, t1.v3.y, t1.v3.z}},
                           K::Triangle_3{K::Point_3{t2.v1.x, t2.v1.y, t2.v1.z},K::Point_3{t2.v2.x, t2.v2.y, t2.v2.z},
                                         K::Point_3{t2.v3.x, t2.v3.y, t2.v3.z}});
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
        /**
        TODO
        */
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
            // couper toutes les edges avec un vertex
            for(Edge e : meshIco->edges()){
                Vertex v0 = e.firstVertex();
                Vertex v1 = e.secondVertex();
                Vector3 vPos0 = geometryIco->inputVertexPositions[v0];
                Vector3 vPos1 = geometryIco->inputVertexPositions[v1];

                Halfedge he = meshIco->insertVertexAlongEdge(e);

                geometryIco->inputVertexPositions[he.vertex()] = (vPos0 + vPos1) * 0.5f / norm((vPos0 + vPos1) * 0.5f);
            }

            //iterer sur les faces pour reconstruire les faces
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

        return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>(std::move(meshIco), std::move(geometryIco));
    }

    Eigen::Matrix3d eulerAngles(double psi, double phi, double theta) {
        return Eigen::Matrix3d{{cos(psi)*cos(phi)-sin(psi)*cos(theta)*sin(phi), -cos(psi)*sin(phi)-sin(psi)*cos(theta)*cos(phi), sin(psi)*sin(theta)},
                               {sin(psi)*cos(phi)+cos(psi)*cos(theta)*sin(phi), -sin(psi)*sin(phi)+cos(psi)*cos(theta)*cos(phi), -cos(psi)*sin(theta)},
                               {sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta)}};
    }

    std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
        mixMeshes(const VertexPositionGeometry& g1, const VertexPositionGeometry& g2) {

        std::vector<Vector3> pts;

        for (const Vertex& v1 : g1.mesh.vertices()) {
            pts.push_back(g1.vertexPositions[v1]);
        }

        std::unordered_map<Vertex, size_t> g2Tog1;
        std::set<Vertex> shared;
        double threshold = 10e-6;
        for (const Vertex& v2 : g2.mesh.vertices()) {
            bool ok = true;
            for (const Vertex& v1 : g1.mesh.vertices()) {
                if (norm(g1.vertexPositions[v1] - g2.vertexPositions[v2]) < threshold) {
                    g2Tog1[v2] = v1.getIndex();
                    shared.insert(v2);
                    ok = false;
                    break;
                }
            }
            if (ok) {
                g2Tog1[v2] = pts.size();
                pts.push_back(g2.vertexPositions[v2]);
            }
        }

        std::vector<std::vector<size_t>> faces;
        for (const Face& f : g1.mesh.faces()) {
            faces.push_back({f.halfedge().vertex().getIndex(), f.halfedge().next().vertex().getIndex(), f.halfedge().next().next().vertex().getIndex()});
        }

        std::set<Vertex> boundary;
        for (const Vertex& v : g2.mesh.boundaryLoop(0).adjacentVertices()) {
            boundary.insert(v);
        }

        /*for (const Face& f : g2.mesh.faces()) {

            if (boundary.count(f.halfedge().vertex()) > 0) continue;
            if (boundary.count(f.halfedge().next().vertex()) > 0) continue;
            if (boundary.count(f.halfedge().next().next().vertex()) > 0) continue;

            faces.push_back({g2Tog1[f.halfedge().vertex()], g2Tog1[f.halfedge().next().vertex()],
                             g2Tog1[f.halfedge().next().next().vertex()]});
        }*/

        for (const Face& f : g2.mesh.faces()) {
            /*int count = 0;
            if (shared.count(f.halfedge().vertex()) > 0) count++;
            if (shared.count(f.halfedge().next().vertex()) > 0) count++;
            if (shared.count(f.halfedge().next().next().vertex()) > 0) count++;
            if (boundary) continue;*/
            if (boundary.count(f.halfedge().vertex()) > 0) continue;
            if (boundary.count(f.halfedge().next().vertex()) > 0) continue;
            if (boundary.count(f.halfedge().next().next().vertex()) > 0) continue;

            faces.push_back({g2Tog1[f.halfedge().vertex()], g2Tog1[f.halfedge().next().vertex()],
            g2Tog1[f.halfedge().next().next().vertex()]});
        }


        return makeManifoldSurfaceMeshAndGeometry(faces, pts);
    }

    void initPCA(VertexPositionGeometry& geom, VertexPositionGeometry& geomTarget, const MeshInfo& meshData, const std::unordered_map<size_t, size_t>& matching) {

        std::vector<Vector3> data;
        for (const Vertex& v : geomTarget.mesh.boundaryLoop(0).adjacentVertices()) {
            data.push_back(geomTarget.vertexPositions[v]);
        }
        Plane p = PCA(data);

        Vector3 centerTarget = Vector3::zero();
        for (const Vector3& v : data)
            centerTarget += v;
        centerTarget /= data.size();

        Vector3 centerSource = Vector3::zero();
        for (const Vertex& v : geom.mesh.vertices()) {
            Vector3 v3 = geom.vertexPositions[v];
            centerSource += v3;
        }
        centerSource /= geom.mesh.nVertices();

        Eigen::Matrix3d R = findRotationToFitVectors(Eigen::Vector3d{0.0, 1.0, 0.0}, Eigen::Vector3d{p.n.x, p.n.y, p.n.z});

        for (const Vertex& v : geom.mesh.vertices()) {
            Eigen::Vector3d vEigen = Eigen::Vector3d{geom.vertexPositions[v].x, geom.vertexPositions[v].y, geom.vertexPositions[v].z};
            vEigen = R * vEigen;
            geom.vertexPositions[v] = Vector3{vEigen[0], vEigen[1], vEigen[2]} - centerSource;
        }

        double theta = optimizeAngle(geom, geomTarget, p.n, centerTarget - centerSource, meshData, matching);
        Eigen::Matrix3d Rtheta;
        Rtheta << cos(theta)+p.n.x*p.n.x*(1-cos(theta)),       p.n.x*p.n.y*(1-cos(theta))-p.n.z*sin(theta), p.n.x*p.n.z*(1-cos(theta))+p.n.y*sin(theta),
                  p.n.y*p.n.x*(1-cos(theta))+p.n.z*sin(theta), cos(theta)+p.n.y*p.n.y*(1-cos(theta)),       p.n.y*p.n.z*(1-cos(theta))-p.n.x*sin(theta),
                  p.n.z*p.n.x*(1-cos(theta))-p.n.y*sin(theta), p.n.z*p.n.y*(1-cos(theta))+p.n.x*sin(theta), cos(theta)+p.n.z*p.n.z*(1-cos(theta));
        for (const Vertex& v : geom.mesh.vertices()) {
            Eigen::Vector3d vEigen = Eigen::Vector3d{geom.vertexPositions[v].x, geom.vertexPositions[v].y, geom.vertexPositions[v].z};
            vEigen = Rtheta * vEigen;
            geom.vertexPositions[v] = 0.5*Vector3{vEigen[0], vEigen[1], vEigen[2]} + centerTarget;
        }

        for (const Vertex& v : meshData.border) {
            size_t id3d = matching.at(v.getIndex());
            geom.vertexPositions[v] = geomTarget.vertexPositions[id3d];
        }

        polyscope::registerSurfaceMesh("PCA", geom.vertexPositions, geom.mesh.getFaceVertexList());
    }

    Plane PCA(std::vector<Vector3>& data) {
        Eigen::Map<Eigen::Matrix3Xd> mat(&data[0][0], 3, data.size());
        Eigen::MatrixX3d matT = mat.transpose();

        Eigen::MatrixX3d dataBar = matT.rowwise() - matT.colwise().mean();
        Eigen::MatrixX3d cov = (dataBar.adjoint() * dataBar) / double(mat.rows() - 1);

        // Find the smallest eigenvalue and associated vector
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
        Eigen::Vector3d eigenValues = es.eigenvalues();
        double min = std::numeric_limits<double>::max(); int id = 0;
        for(int j = 0; j < 3; ++j){
            if(eigenValues[j] < min){
                min = eigenValues[j];
                id = j;
            }
        }
        Vector3 n = Vector3{es.eigenvectors().col(id)[0],
                            es.eigenvectors().col(id)[1],
                            es.eigenvectors().col(id)[2]};
        n = n.normalize();

        return Plane(-n);
    }

    Eigen::Matrix3d findRotationToFitVectors(Eigen::Vector3d a, Eigen::Vector3d b){
        a = a.normalized();
        b = b.normalized();

        Eigen::Vector3d v = a.cross(b);
        double c = a.dot(b);
        Eigen::Matrix3d Vcross;
        Vcross << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + Vcross + Vcross*Vcross*(1.0/(1.0 + c));
        return R;
    }

    double optimizeAngle(VertexPositionGeometry& geom, VertexPositionGeometry& geomTarget, const Vector3& u, const Vector3& shift,
                         const MeshInfo& meshData, const std::unordered_map<size_t, size_t>& matching) {
        constexpr int N = 128;
        double bestTheta = 0.0;
        double bestE = std::numeric_limits<double>::max();

        for (int i = 0; i < N; ++i) {
            double theta = i * 2*M_PI/N;
            double E = 0.0;
            Eigen::Matrix3d R;
            R << cos(theta)+u.x*u.x*(1-cos(theta)),     u.x*u.y*(1-cos(theta))-u.z*sin(theta), u.x*u.z*(1-cos(theta))+u.y*sin(theta),
                 u.y*u.x*(1-cos(theta))+u.z*sin(theta), cos(theta)+u.y*u.y*(1-cos(theta)),     u.y*u.z*(1-cos(theta))-u.x*sin(theta),
                 u.z*u.x*(1-cos(theta))-u.y*sin(theta), u.z*u.y*(1-cos(theta))+u.x*sin(theta), cos(theta)+u.z*u.z*(1-cos(theta));
            for (const Vertex& v : meshData.border) {
                Eigen::Vector3d vEigen = Eigen::Vector3d{geom.vertexPositions[v].x, geom.vertexPositions[v].y, geom.vertexPositions[v].z};
                vEigen = R * vEigen;
                Vector3 vec = Vector3{vEigen[0], vEigen[1], vEigen[2]} + shift;
                size_t id3d = matching.at(v.getIndex());
                E = norm(geomTarget.vertexPositions[id3d] - vec);
            }

            if (E < bestE) {
                bestTheta = theta;
                bestE = E;
            }
        }
        bestTheta = 3.18;
        std::cout << "bestTheta : " << bestTheta << std::endl;
        return bestTheta;
    }

    Eigen::Vector3d to_eigen(const Vector3& v) {
        return Eigen::Vector3d{v.x, v.y, v.z};
    }

    Vector3 to_geometrycentral(const Eigen::Vector3d& v) {
        return Vector3{v[0], v[1], v[2]};
    }

}