#include "application.hpp"

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
LBFGSopt* Application::lbfgs = nullptr;
TinyADopt* Application::tinyOpt = nullptr;
int Application::highlight = 0;
ARAP* Application::arapHolder = nullptr;
int Application::index = 0;
double Application::x = 0;
double Application::y = 0;
double Application::z = 0;

VertexData<std::array<Vector3,2>>* Application::planeVertexTangent = nullptr;

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
            if (S_ISDIR(st.st_mode)) {
                if (ent->d_name[0] != '.')
                    directories.push_back(ent->d_name);
            }
            else if (S_ISREG(st.st_mode)) {
                files.emplace_back(ent->d_name);
            }
        }
        closedir(dir);
    }
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
    std::cout << "#####" << std::endl;
    for (const Halfedge& he : BL2.adjacentHalfedges()) {
        double edge = norm(geometry2D->vertexPositions[he.next().vertex()] - geometry2D->vertexPositions[he.vertex()]);
        double edge0 = edgeLengths[matching[he.vertex().getIndex()]];
        double angle = geometry2D->vertexAngleSums[he.vertex()];
        double angle0 = angles[matching[he.vertex().getIndex()]];

        Eedge += (edge - edge0)*(edge - edge0);
        Eangle += (angle - angle0)*(angle - angle0);
        N++;
        std::cout << "edge0 " << matching[he.vertex().getIndex()] << " : " << edge0 << std::endl;
        std::cout << "edge " << he.vertex() << " : " << edge << std::endl;
        std::cout << "angle0 " << matching[he.vertex().getIndex()] << " : " << angle0 << std::endl;
        std::cout << "angle " << he.vertex() << " : " << angle << std::endl;
    }
    geometry2D->unrequireVertexAngleSums();

    std::cout << "Eangle : " << Eangle/N << std::endl;
    std::cout << "Eedge : " << Eedge/N << std::endl;
}

void Application::make2DMesh() {

    //make2DLarge();
    //return;

    BoundaryLoop BL = mesh->boundaryLoop(0);//2

    std::map<size_t, double> angles;
    std::map<size_t, double> edgeLengths;
    geometry->requireVertexAngleSums();
    for (const Halfedge& he : BL.adjacentHalfedges()) {
        edgeLengths[he.vertex().getIndex()] = norm(geometry->vertexPositions[he.next().vertex()] - geometry->vertexPositions[he.vertex()]);
        angles[he.vertex().getIndex()] = geometry->vertexAngleSums[he.vertex()];
        std::cout << "edge" << he.vertex().getIndex() << " : " << edgeLengths[he.vertex().getIndex()] << std::endl;
        std::cout << "angle " << he.vertex().getIndex() << " : " << angles[he.vertex().getIndex()] << std::endl;
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

    for (std::vector<size_t> elem : faces) {
        std::cout << elem[0] << ", " << elem[1] << ", " << elem[2] << std::endl;
        std::cout << matching[elem[0]] << ", " << matching[elem[1]] << ", " << matching[elem[2]] << std::endl;
    }

    //faces[faces.size()-1] = {{1, 0, 486}};

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
    bc<<0,0,l,0;

    // LSCM parametrization
    igl::lscm(V,F,b,bc,V_uv);

    //ARAP Param
    /*Eigen::MatrixXd initial_guess;
    Eigen::VectorXi bnd;
    igl::boundary_loop(F,bnd);
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V,bnd,bnd_uv);

    igl::harmonic(V,F,bnd,bnd_uv,1,initial_guess);

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 1000;
    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V,F,2,b,arap_data);


    // Solve arap using the harmonic map as initial guess
    V_uv = initial_guess;

    arap_solve(bc,arap_data,V_uv);*/

    Eigen::MatrixXd V_uvdef(V.rows(),3);
    for (Eigen::Index i = 0; i < V_uv.rows(); ++i) {
        V_uvdef.coeffRef(i, 0) = V_uv.coeffRef(i,0);
        V_uvdef.coeffRef(i, 1) = 0;
        V_uvdef.coeffRef(i, 2) = V_uv.coeffRef(i,1);
    }

    for (Eigen::Index i = 0; i < V_uv.rows(); ++i) {
        pts[i] = Vector3{V_uv.coeffRef(i,0), 0, V_uv.coeffRef(i,1)};
    }
    //polyscope::registerSurfaceMesh("param mesh", V_uvdef, F);




    std::tie(mesh2D, geometry2D) = makeManifoldSurfaceMeshAndGeometry(faces, pts);
    show2DMetrics(edgeLengths, angles);
    writeSurfaceMesh(*mesh2D, *geometry2D, "../../data/2Dmesh.obj");
    std::cout << "mesh wrote" << std::endl;
    //polyscope::registerSurfaceMesh("Parametrisation", geometry2D->vertexPositions, mesh2D->getFaceVertexList());
}

void Application::init() {
    polyscope::init();
    polyscope::state::userCallback = callbacks[0];
    TraverseDirectory("../../data");

    polyscope::show();
}

void Application::readOBJ(const std::string& fn) {
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(fn);

    for (const Vertex& v : mesh->vertices()) {
        //geometry->vertexPositions[v] += 0.015*Vector3{polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()};
    }

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


        std::vector<std::array<double, 3>> fColor(mesh->nFaces());
        for (size_t iF = 0; iF < mesh->nFaces(); iF++) {
            Face f = mesh->face(iF);
            bool boundary = false;
            if (f.halfedge().vertex().isBoundary())
                boundary = true;
            if (f.halfedge().next().vertex().isBoundary())
                boundary = true;
            if (f.halfedge().next().next().vertex().isBoundary())
                boundary = true;

            if (!boundary)
                fColor[iF] = {{28.0/255, 99.0/255, 227.0/255}};
            else
                fColor[iF] = {{99.0/255, 227.0/255, 28.0/255}};
        }
        polyscope::getSurfaceMesh("Mesh")->addFaceColorQuantity("fColor", fColor);
        polyscope::registerSurfaceMesh("Border", geometry2D->vertexPositions, mesh2D->getFaceVertexList());


        //polyscope::registerSurfaceMesh("Parametrisation", geometry2D->vertexPositions, mesh2D->getFaceVertexList());
        triangulateLiepa(geometry2D, mesh2D, Utils::averageBlEdgeLength(mesh2D->boundaryLoop(1), *geometry2D));
        //triangulateLiepaFair(geometry, mesh);

        meshData2D = new MeshInfo(*mesh2D, *geometry2D);
        meshData2D->precomputeData();
        geometry2D->requireVertexTangentBasis();
        planeVertexTangent = new VertexData<std::array<Vector3,2>>(*mesh2D);
                geometry2D->vertexTangentBasis;
        for (Vertex v : mesh2D->vertices()) {
            (*planeVertexTangent)[v] = std::array<Vector3,2>{geometry2D->vertexTangentBasis[v][0], geometry2D->vertexTangentBasis[v][1]};
        }

        geometry2D->unrequireVertexTangentBasis();
        polyscope::registerSurfaceMesh("Parametrisation", geometry2D->vertexPositions, mesh2D->getFaceVertexList());

        std::vector<std::array<double, 3>> fColor2(mesh2D->nFaces());
        for (size_t iF = 0; iF < mesh2D->nFaces(); iF++) {
            Face f = mesh2D->face(iF);
            bool boundary = false;
            if (f.halfedge().vertex().isBoundary())
                boundary = true;
            if (f.halfedge().next().vertex().isBoundary())
                boundary = true;
            if (f.halfedge().next().next().vertex().isBoundary())
                boundary = true;

            if (!boundary)
                fColor2[iF] = {{28.0/255, 99.0/255, 227.0/255}};
            else
                fColor2[iF] = {{99.0/255, 227.0/255, 28.0/255}};
        }
        polyscope::getSurfaceMesh("Parametrisation")->addFaceColorQuantity("fColor2", fColor2);

        //ARAP arapHolder(*geometry2D);
        //Utils::harmonicSurface(*geometry2D, *geometry, *meshData2D, matching);
        arapHolder = new ARAP(*geometry2D);
        Utils::ARAP(*geometry2D, *geometry, *meshData2D, matching, 0.5);
        //Utils::initPCA(*geometry2D, *geometry, *meshData2D, matching);
        gd = new GradientDescent(*meshData2D);
        lbfgs = new LBFGSopt(*meshData2D);
        tinyOpt = new TinyADopt(*meshData2D);

        polyscope::registerSurfaceMesh("Optim", geometry2D->vertexPositions, mesh2D->getFaceVertexList());
        polyscope::getSurfaceMesh("Optim")->addVertexScalarQuantity("radius", meshData2D->repulsiveRadius);

        std::unique_ptr<ManifoldSurfaceMesh> m;
        std::unique_ptr<VertexPositionGeometry> g;
        //std::tie(m, g) = readManifoldSurfaceMesh("../../data/caryatid.obj");
        //polyscope::registerSurfaceMesh("Bonus mesh", g->vertexPositions, m->getFaceVertexList());

        changeState(OPTI_1D);
    }
}

void Application::callback1() {

    ImGui::SliderFloat("Alpha : ", &gd->alpha, 0.0f, 1.0f);
    ImGui::SliderFloat("Beta : ", &gd->beta, 0.0f, 0.1f);
    ImGui::SliderFloat("Gamma : ", &gd->gamma, 0.0f, 0.1f);
    ImGui::SliderFloat("Lambda : ", &gd->lambda, 0.0f, 1.0f);
    ImGui::SliderFloat("Mu : ", &gd->mu, 0.0f, 1.0f);
    ImGui::SliderFloat("Time step : ", &gd->timeStep, 0.0f, 1.0f);

    if (ImGui::Button("Step")) {
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        for (int j = 0; j <= 100; ++j) {
            //gd->lambda = 0.5f + j*0.005;
            for (int i = 0; i < 50; ++i) {
                gd->step();
            }
        }
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "step() : " << elapsed_seconds.count() << std::endl;
    }
    if (ImGui::Button("Toggle")) {
        toggle = !toggle;
    }
    if (ImGui::Button("Show")) {
        polyscope::getSurfaceMesh("Optim")->updateVertexPositions(meshData2D->geom.vertexPositions);
    }
    if (ImGui::Button("LBFGS")) {
        lbfgs->run();
    }
    if (ImGui::Button("TinyAD")) {
        tinyOpt->run();
    }
    if (toggle) {
        /*if (x != 0.0 )
            meshData2D->geom.vertexPositions[index] = Vector3{x,y,z};*/
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        gd->step();
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "step() : " << elapsed_seconds.count() << std::endl;
    }
    ImGui::Separator();
    ImGui::InputInt("Vertex : ", &highlight);
    if (ImGui::Button("Highlight : ")) {
        VertexData<bool> highlightData(*mesh2D, false);
        highlightData[highlight] = true;
        polyscope::getSurfaceMesh("ARAP")->addVertexScalarQuantity("highlight", highlightData);
    }
    ImGui::Separator();
    if (ImGui::Button("ARAP")) {
        Utils::ARAP(*arapHolder, *geometry2D, *geometry, *meshData2D, matching);
    }
    if (ImGui::Button("Mix")) {
        std::unique_ptr<ManifoldSurfaceMesh> m;
        std::unique_ptr<VertexPositionGeometry> g;
        std::tie(m, g) = Utils::mixMeshes(*geometry, *geometry2D);
        writeSurfaceMesh(*m, *g, "../../data/surface.obj");
        polyscope::registerSurfaceMesh("Mix", g->vertexPositions, m->getFaceVertexList());
    }
    if (ImGui::Button("Write")) {
        writeSurfaceMesh(meshData2D->mesh, meshData2D->geom, "../../data/surface.obj");
    }
    if (ImGui::Button("Ruling")) {

        geometry2D->requireVertexTangentBasis();

        VertexData<curvatureData> data = curvatureProcessing::calculateCurvature(*geometry2D, 0.1);
        VertexData<Vector3> u1(*mesh2D);
        for (const Vertex& v : mesh2D->vertices()) {
            u1[v] = data[v].u1;
            if (v.isBoundary())
                u1[v] = Vector3::zero();
        }

        VertexData<Vector3> u1_proj(*mesh2D);
        VertexData<Vector2> u1_2D(*mesh2D);
        VertexData<Vector3> vBasisX(*mesh2D);
        for(Vertex v : mesh2D->vertices()){
            Eigen::Vector3d u; u << u1[v].x, u1[v].y, u1[v].z;
            Eigen::MatrixXd A(3, 2);
            A << geometry2D->vertexTangentBasis[v][0][0], geometry2D->vertexTangentBasis[v][1][0],
                    geometry2D->vertexTangentBasis[v][0][1], geometry2D->vertexTangentBasis[v][1][1],
                    geometry2D->vertexTangentBasis[v][0][2], geometry2D->vertexTangentBasis[v][1][2];

            auto b = (A.transpose() * A).inverse() * A.transpose() * u;
            Vector2 u1_proj_2D = Vector2{b[0], b[1]};
            u1_2D[v] = u1_proj_2D;
            vBasisX[v] = (*planeVertexTangent)[v][0];
            /*u1_proj[v] = Vector3{ (*planeVertexTangent)[v][0][0] * u1_proj_2D[0] + (*planeVertexTangent)[v][1][0] * u1_proj_2D[1],
                                  (*planeVertexTangent)[v][0][1] * u1_proj_2D[0] + (*planeVertexTangent)[v][1][1] * u1_proj_2D[1],
                                  (*planeVertexTangent)[v][0][2] * u1_proj_2D[0] + (*planeVertexTangent)[v][1][2] * u1_proj_2D[1]};*/
        }

        geometry2D->unrequireVertexTangentBasis();

        polyscope::getSurfaceMesh("Parametrisation")->setVertexTangentBasisX(vBasisX);
        polyscope::getSurfaceMesh("Parametrisation")->addVertexIntrinsicVectorQuantity("ruling dirs", u1_2D);
        polyscope::getSurfaceMesh("Optim")->addVertexVectorQuantity("ruling dirs", u1);
    }
}




