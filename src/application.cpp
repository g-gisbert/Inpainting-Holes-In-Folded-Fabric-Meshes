#include "application.hpp"

#include <utility>

State Application::state = CHOOSE_FILE;
std::unique_ptr<ManifoldSurfaceMesh> Application::mesh;
std::unique_ptr<VertexPositionGeometry> Application::geometry;
std::unique_ptr<ManifoldSurfaceMesh> Application::mesh2D;
std::unique_ptr<VertexPositionGeometry> Application::geometry2D;
void (*Application::callbacks[2])() = { Application::callback0, Application::callback1 };
std::string Application::currentPath, Application::filename;
std::vector<std::string> Application::files, Application::directories;
MeshInfo* Application::meshData2D = nullptr;
std::unordered_map<size_t, size_t> Application::matching;
bool Application::toggle = false;
GradientDescent* Application::gd = nullptr;
int Application::highlight = 0;
int Application::index = 0;
bool Application::bendingEnergy = true;
bool Application::radioOption = false;
bool Application::enableCollisions = false;


void Application::TraverseDirectory(const std::string& path) {
    currentPath = path;
    files.clear();
    directories.clear();
    directories.emplace_back("..");
    DIR* dir = opendir(path.c_str());
    if (dir) {
        struct dirent* ent;
        while ((ent = readdir(dir)) != nullptr) {
            struct stat st{};
            std::string fullpath = path + "/" + ent->d_name;
            if (stat(fullpath.c_str(), &st) == -1)
                continue;
            if (S_ISREG(st.st_mode)) {
                files.emplace_back(ent->d_name);
            }
        }
        closedir(dir);
    }
    std::sort(files.begin(), files.end());
}

void Application::changeState(State newState) {
    polyscope::state::userCallback = callbacks[newState];
}

void Application::show2DMetrics(std::map<size_t, double>& edgeLengths, std::map<size_t, double>& angles) {
    geometry2D->requireVertexAngleSums();
    BoundaryLoop BL2 = mesh2D->boundaryLoop(1);
    double Eangle = 0.0;
    double Eedge = 0.0;
    int N = 0;
    std::cout << "2D parametrisation" << std::endl;
    for (const Halfedge& he : BL2.adjacentHalfedges()) {
        double edge = norm(geometry2D->vertexPositions[he.next().vertex()] - geometry2D->vertexPositions[he.vertex()]);
        double edge0 = edgeLengths[matching[he.vertex().getIndex()]];
        double angle = geometry2D->vertexAngleSums[he.vertex()];
        double angle0 = angles[matching[he.vertex().getIndex()]];

        Eedge += (edge - edge0)*(edge - edge0);
        Eangle += (angle - angle0)*(angle - angle0);
        N++;
    }
    geometry2D->unrequireVertexAngleSums();

    std::cout << "Eangle : " << Eangle/N << std::endl;
    std::cout << "Eedge : " << Eedge/N << std::endl;
}

void Application::make2DMesh() {

    BoundaryLoop BL = mesh->boundaryLoop(0);

    std::map<size_t, double> angles;
    std::map<size_t, double> edgeLengths;
    geometry->requireVertexAngleSums();
    for (const Halfedge& he : BL.adjacentHalfedges()) {
        edgeLengths[he.vertex().getIndex()] = norm(geometry->vertexPositions[he.next().vertex()] - geometry->vertexPositions[he.vertex()]);
        angles[he.vertex().getIndex()] = geometry->vertexAngleSums[he.vertex()];
    }
    geometry->unrequireVertexAngleSums();

    std::unordered_set<size_t> BLvertices;
    for (Vertex v : BL.adjacentVertices()) {
        BLvertices.insert(v.getIndex());
    }

    std::vector<Vector3> pts;
    std::vector<std::vector<size_t>> faces;

    Halfedge startHe = BL.halfedge().next().twin().next();
    Halfedge he = startHe;
    pts.push_back(geometry->vertexPositions[startHe.vertex()]);
    pts.push_back(geometry->vertexPositions[startHe.tipVertex()]);
    matching[0] = startHe.vertex().getIndex();
    matching[1] = startHe.tipVertex().getIndex();
    size_t id0 = startHe.vertex().getIndex();
    size_t id1 = startHe.tipVertex().getIndex();
    size_t lastIdBoundary = 0;
    size_t lastIdIn = 1;
    size_t vertexNumber = 1;

    do {
        he = he.next();
        Vertex v = he.tipVertex();

        if (v.getIndex() == id0)
            faces.push_back({lastIdBoundary, lastIdIn, 0});
        else if (v.getIndex() == id1)
            faces.push_back({lastIdBoundary, lastIdIn, 1});
        else {
            pts.push_back(geometry->vertexPositions[v]);
            matching[++vertexNumber] = v.getIndex();
            faces.push_back({lastIdBoundary, lastIdIn, vertexNumber});
        }
        if (BLvertices.count(v.getIndex()) == 1) {
            lastIdBoundary = vertexNumber;
            if (v.getIndex() == id0)
                lastIdBoundary = 0;
            if (v.getIndex() == id1)
                lastIdBoundary = 1;
            he = he.twin();
        } else {
            lastIdIn = vertexNumber;
            if (v.getIndex() == id0)
                lastIdBoundary = 0;
            if (v.getIndex() == id1)
                lastIdBoundary = 1;
            he = he.next().twin();
        }
    } while(startHe != he);


    std::tie(mesh2D, geometry2D) = makeManifoldSurfaceMeshAndGeometry(faces, pts);

    Eigen::MatrixXd V(pts.size(), 3);
    Eigen::MatrixXi F(faces.size(), 3);
    Eigen::MatrixXd V_uv;

    for (size_t i = 0; i < pts.size(); ++i) {
        V(i, 0) = pts[i][0];
        V(i, 1) = pts[i][1];
        V(i, 2) = pts[i][2];
    }
    for (size_t i = 0; i < faces.size(); ++i) {
        F(i, 0) = faces[i][0];
        F(i, 1) = faces[i][1];
        F(i, 2) = faces[i][2];
    }

    //LSCM
    Eigen::MatrixXi b(2,1);
    b(0) = 0;
    b(1) = 1;
    double l = norm(pts[0] - pts[1]);
    Eigen::MatrixXd bc(2,2);
    bc << 0, 0, l, 0;

    // LSCM parametrization
    igl::lscm(V,F,b,bc,V_uv);

    Eigen::MatrixXd V_uvdef(V.rows(),3);
    for (Eigen::Index i = 0; i < V_uv.rows(); ++i) {
        V_uvdef.coeffRef(i, 0) = V_uv.coeffRef(i,0);
        V_uvdef.coeffRef(i, 1) = 0;
        V_uvdef.coeffRef(i, 2) = V_uv.coeffRef(i,1);
    }

    for (Eigen::Index i = 0; i < V_uv.rows(); ++i) {
        pts[i] = Vector3{V_uv.coeffRef(i,0), 0, V_uv.coeffRef(i,1)};
    }

    std::tie(mesh2D, geometry2D) = makeManifoldSurfaceMeshAndGeometry(faces, pts);
    show2DMetrics(edgeLengths, angles);
}

void Application::init(const std::string& df) {
    polyscope::init();
    polyscope::state::userCallback = callbacks[0];
    TraverseDirectory(df);

    polyscope::show();
}

void Application::readOBJ(const std::string& fn) {
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(fn);
}

void Application::callback0() {
    ImGui::Text("Current path: %s", currentPath.c_str());
    ImGui::Separator();
    ImGui::Text("Directories:");
    for (const std::string& dir : directories) {
        if (ImGui::Selectable(dir.c_str())) {
            TraverseDirectory(currentPath + "/" + dir);
        }
    }
    ImGui::Separator();
    ImGui::Text("Files:");
    for (const std::string& file : files) {
        if (ImGui::Selectable(file.c_str())) {
            filename = file;
        }
    }
    ImGui::Separator();
    if (ImGui::Button("Choose")) {
        readOBJ(currentPath + "/" + filename);
        make2DMesh();
        polyscope::registerSurfaceMesh("Mesh", geometry->vertexPositions, mesh->getFaceVertexList());

        polyscope::registerSurfaceMesh("Border", geometry2D->vertexPositions, mesh2D->getFaceVertexList());

        triangulate(geometry2D, mesh2D, Utils::averageBlEdgeLength(mesh2D->boundaryLoop(1), *geometry2D));

        meshData2D = new MeshInfo(*mesh2D, *geometry2D);
        meshData2D->precomputeData();
        polyscope::registerSurfaceMesh("Parametrisation", geometry2D->vertexPositions, mesh2D->getFaceVertexList());

        Utils::ARAP(*geometry2D, *geometry, *meshData2D, matching, 0.75);
        gd = new GradientDescent(*meshData2D);

        polyscope::registerSurfaceMesh("Optim", geometry2D->vertexPositions, mesh2D->getFaceVertexList());
        polyscope::getSurfaceMesh("Optim")->addVertexScalarQuantity("radius", meshData2D->repulsiveRadius);

        changeState(OPTI_1D);
    }
}

void Application::callback1() {

    ImGui::SliderFloat("Edge Lengths : ", &gd->alpha, 0.0f, 1.0f);
    ImGui::SliderFloat("Bending : ", &gd->beta, 0.0f, 0.1f);
    ImGui::SliderFloat("Repulsive energy : ", &gd->mu, 0.0f, 1.0f);
    ImGui::SliderFloat("Time step : ", &gd->timeStep, 0.0f, 1.0f);

    if (ImGui::RadioButton("Quadratic Bending Energy", bendingEnergy))
    {
        bendingEnergy = true;
        radioOption = false;
    }

    if (ImGui::RadioButton("2-Ring Edge Lengths", radioOption))
    {
        bendingEnergy = false;
        radioOption = true;
    }
    if (ImGui::Checkbox("Enable collisions", &enableCollisions)) {
        std::cout << "enableCollisions : " << enableCollisions << std::endl;
    }

    ImGui::Separator();
    if (ImGui::Button("Steps")) {
        for (int j = 0; j <= 1000; ++j) {
            gd->step(bendingEnergy, enableCollisions);
        }
        std::cout << "ok";
        polyscope::getSurfaceMesh("Optim")->addVertexScalarQuantity("abc", meshData2D->repulsiveRadius);
        polyscope::getSurfaceMesh("Optim")->updateVertexPositions(meshData2D->geom.vertexPositions);
    }
    if (ImGui::Button("Toggle")) {
        toggle = !toggle;
    }
    if (ImGui::Button("Update Mesh")) {
        //polyscope::getSurfaceMesh("Optim")->updateVertexPositions(meshData2D->geom.vertexPositions);
    }

    if (toggle) {
        gd->step(bendingEnergy, enableCollisions);
    }
    ImGui::Separator();
    if (ImGui::Button("Save mesh")) {
        writeSurfaceMesh(meshData2D->mesh, meshData2D->geom, "surface.obj");
    }

}




