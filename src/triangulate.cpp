#include "triangulate.hpp"


void triangulate(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh, double size) {
    CDT cdt;

    std::vector<Vertex_handle> vertices;
    BoundaryLoop BL = mesh->boundaryLoop(1);
    for(Vertex v : BL.adjacentVertices()){
        Vector3 p = geom->vertexPositions[v];
        Vertex_handle vh = cdt.insert(Point(p.x, p.z));
        vertices.push_back(vh);
    }
    for(size_t i = 0; i < vertices.size(); ++i){
        cdt.insert_constraint(vertices[i], vertices[(i+1)%int(vertices.size())]);
    }

    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, size));
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

    std::map<Vertex_handle, std::size_t> verticesIndex;
    size_t index = 0;
    for (CDT::Finite_vertices_iterator t = cdt.vertices_begin(); t != cdt.vertices_end(); t++) {
        verticesIndex[t] = index++;
    }

    std::unique_ptr<ManifoldSurfaceMesh> meshTmp;
    std::unique_ptr<VertexPositionGeometry> geometryTmp;
    std::tie(meshTmp, geometryTmp) = toGC(cdt, verticesIndex);
    std::tie(meshTmp, geometryTmp) = cleanMesh(*geom, *geometryTmp);

    std::tie(meshTmp, geometryTmp) = mixMeshes(*geom, *geometryTmp);

    mesh = std::move(meshTmp);
    geom = std::move(geometryTmp);
}

void triangulateLiepa(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh, double size) {

    writeSurfaceMesh(*mesh, *geom, "tmp.obj");

    Mesh meshCGAL;
    if(!PMP::IO::read_polygon_mesh(std::string("tmp.obj"), meshCGAL))
    {
        std::cerr << "Invalid input." << std::endl;
        return;
    }

    std::vector<halfedge_descriptor> border_cycles;
    // collect one halfedge per boundary cycle
    PMP::extract_boundary_cycles(meshCGAL, std::back_inserter(border_cycles));
    halfedge_descriptor h = border_cycles[0];
    std::vector<face_descriptor>  patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    PMP::triangulate_and_refine_hole(meshCGAL,
        h, std::back_inserter(patch_facets), std::back_inserter(patch_vertices),
        CGAL::parameters::density_control_factor(1.7)); //1.4
    CGAL::IO::write_polygon_mesh("tmp2.obj", meshCGAL, CGAL::parameters::stream_precision(17));
    std::tie(mesh, geom) = readManifoldSurfaceMesh("tmp2.obj");
}


void triangulateLiepaFair(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh) {

    writeSurfaceMesh(*mesh, *geom, "tmp.obj");
    size_t N = mesh->nVertices();

    Mesh meshCGAL;
    if(!PMP::IO::read_polygon_mesh(std::string("tmp.obj"), meshCGAL))
    {
        std::cerr << "Invalid input." << std::endl;
        return;
    }

    std::vector<halfedge_descriptor> border_cycles;
    // collect one halfedge per boundary cycle
    PMP::extract_boundary_cycles(meshCGAL, std::back_inserter(border_cycles));
    halfedge_descriptor h = border_cycles[0];
    std::vector<face_descriptor>  patch_facets;
    std::vector<vertex_descriptor> patch_vertices;

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    PMP::triangulate_refine_and_fair_hole(meshCGAL,
                                          h, std::back_inserter(patch_facets), std::back_inserter(patch_vertices),
                                          CGAL::parameters::density_control_factor(1.7));


    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "liepa() : " << elapsed_seconds.count() << std::endl;

    CGAL::IO::write_polygon_mesh("tmp2.obj", meshCGAL, CGAL::parameters::stream_precision(17));
    std::unique_ptr<ManifoldSurfaceMesh> m;
    std::unique_ptr<VertexPositionGeometry> g;
    std::tie(m, g) = readManifoldSurfaceMesh("tmp2.obj");

    // AT
    std::ofstream file("../../data/file.obj");
    for (size_t i = 0; i < N; ++i) {
        Vertex v = m->vertex(i);
        file << "v " << g->vertexPositions[v].x << " " << g->vertexPositions[v].y << " " << g->vertexPositions[v].z << " 1.000000 1.000000 1.000000\n";
    }
    for (size_t i = N; i < m->nVertices(); ++i) {
        Vertex v = m->vertex(i);
        file << "v " << g->vertexPositions[v].x << " " << g->vertexPositions[v].y << " " << g->vertexPositions[v].z << " 0.000000 1.000000 0.000000\n";
    }

    for (const std::vector<size_t>& f : m->getFaceVertexList()) {
        file << "f " << f[0]+1 << " " << f[1]+1 << " " << f[2]+1 << "\n";
    }
    file.close();


    polyscope::registerSurfaceMesh("Liepa", g->vertexPositions, m->getFaceVertexList());
}


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
toGC(CDT& cdt, std::map<Vertex_handle, std::size_t>& verticesIndex){

    int i = 0;
    std::vector<Vector3> pts(cdt.number_of_vertices());
    for (CDT::Finite_vertices_iterator t = cdt.vertices_begin(); t != cdt.vertices_end(); t++) {
        pts[i++] = Vector3{t->point()[0], 0, t->point()[1]};
    }

    i = 0;
    std::vector<std::vector<size_t>> faces(cdt.number_of_faces(), std::vector<size_t>(3, 0));
    std::cout << "before : " << faces.size() << std::endl;
    for (CDT::Finite_faces_iterator t = cdt.faces_begin(); t != cdt.faces_end(); t++) {
        if(t->is_in_domain())
            faces[i++] = {size_t(verticesIndex[t->vertex(0)]),
                          size_t(verticesIndex[t->vertex(2)]),
                          size_t(verticesIndex[t->vertex(1)])};
    }
    faces.resize(i);
    std::cout << "total : " << cdt.number_of_faces() << std::endl;
    std::cout << "in domain : " << i << std::endl;
    std::cout << "after : " << faces.size() << std::endl;

    return makeManifoldSurfaceMeshAndGeometry(faces, pts);
}


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
cleanMesh(VertexPositionGeometry& geom, VertexPositionGeometry& geomFlat){

    std::vector<int> vertices(geomFlat.mesh.nVertices());
    std::vector<int> faces(geomFlat.mesh.nFaces());
    std::set<size_t> verticesToDelete;
    std::map<size_t, size_t> correspondance;
    int i = 0;
    for(Vertex v : geomFlat.mesh.vertices()){
        vertices[i++] = v.getIndex();
    }
    i = 0;
    for(Face f : geomFlat.mesh.faces()){
        faces[i++] = f.getIndex();
    }
    std::vector<Vector3> positions(geomFlat.mesh.nVertices());
    for(size_t i = 0; i < geomFlat.mesh.nVertices(); ++i)
        positions[i] = geomFlat.vertexPositions[i];
    std::vector<std::vector<size_t>> faceVertexList = geomFlat.mesh.getFaceVertexList();
    std::vector<size_t> facesIdToDelete;

    for(Vertex v : geomFlat.mesh.boundaryLoop(0).adjacentVertices()){
        for(Vertex vb : geom.mesh.boundaryLoop(1).adjacentVertices()){
            if(abs(norm(geom.vertexPositions[vb] - geomFlat.vertexPositions[v])) < 10e-8){
                correspondance[v.getIndex()] = vb.getIndex();
                break;
            }
        }
    }

    Halfedge current = geomFlat.mesh.boundaryLoop(0).halfedge();
    while(correspondance.count(current.vertex().getIndex()) == 0)
        current = current.next();
    size_t start = current.vertex().getIndex();
    do{
        std::cout << "#########" << std::endl;
        std::cout << "prev vertex id : " << current.vertex().getIndex() << std::endl;
        Halfedge prev = current;
        current = current.next();
        std::cout << "curr vertex id : " << current.vertex().getIndex() << std::endl;
        std::vector<size_t> inBetweenVertices;
        while(true){
            // if bad vertex
            size_t id = current.vertex().getIndex();
            if(correspondance.count(id) == 0){
                vertices[id] = -1;
                verticesToDelete.insert(id);
                inBetweenVertices.push_back(id);
                current = current.next();
                std::cout << "\trolling vertex id : " << current.vertex().getIndex() << std::endl;
            }
            else{
                break;
            }
        }
        std::cout << "size : " << inBetweenVertices.size() << std::endl;
        std::vector<size_t> oneRingToConnect;
        for(size_t j = 0; j < inBetweenVertices.size(); ++j){
            size_t vId = inBetweenVertices[j];
            Vertex v = geomFlat.mesh.vertex(vId);
            for(Face f : v.adjacentFaces()){
                if(std::count(facesIdToDelete.begin(), facesIdToDelete.end(), f.getIndex()) == 0){
                    //facesIdToDelete.push_back(f.getIndex());
                }
            }

            // Good order
            Halfedge it = v.halfedge().next();
            while(true){
                if(j == 0 && it.vertex().getIndex() == prev.vertex().getIndex())
                    break;
                if(j != 0 && it.vertex().getIndex() == inBetweenVertices[j-1])
                    break;
                it = it.next().twin().next();
            }

            for(size_t c = 0; c < v.degree(); ++c){
                size_t id = it.vertex().getIndex();
                if(std::count(inBetweenVertices.begin(), inBetweenVertices.end(), id) == 0
                   && std::count(oneRingToConnect.begin(), oneRingToConnect.end(), id) == 0
                   && prev.vertex().getIndex() != id && current.vertex().getIndex() != id
                   && vertices[id] != -1){
                    oneRingToConnect.push_back(id);
                }
                it = it.next().twin().next();
            }
        }
        oneRingToConnect.insert(oneRingToConnect.begin(), 1, prev.vertex().getIndex());
        for(int c = 0; c < int(oneRingToConnect.size())-1; ++c){
            faceVertexList.push_back({oneRingToConnect[c], oneRingToConnect[c+1], current.vertex().getIndex()});
            std::cout << "added : " << oneRingToConnect[c] << ", "
                      << oneRingToConnect[c+1] << ", " << current.vertex().getIndex() << std::endl;
        }

    }while(current.vertex().getIndex() != start);

    std::cout << "To resolve" << std::endl;
    // Resolve
    i = 0;
    for(std::vector<size_t>& face : faceVertexList){
        if(verticesToDelete.count(face[0]) > 0 ||
           verticesToDelete.count(face[1]) > 0 || verticesToDelete.count(face[2]) > 0){
            facesIdToDelete.push_back(i);
        }
        i++;
    }
    std::sort(facesIdToDelete.begin(), facesIdToDelete.end(), std::greater<size_t>());
    for(size_t id : facesIdToDelete){
        std::cout << "deleted face : " << id << std::endl;
        faceVertexList.erase(faceVertexList.begin() + id);
    }
    std::cout << "shifting" << std::endl;
    std::map<size_t, size_t> match;
    std::vector<Vector3> points;
    size_t shift = 0;
    for(int id : vertices){
        // shift tout ca
        if(id == -1){
            shift++;
        }
        else{
            match[id] = id - shift;
            std::cout << id << " -> " << id-shift << std::endl;
            points.push_back(geomFlat.vertexPositions[id]);
        }
    }
    // changer les id de faces
    for(std::vector<size_t>& face : faceVertexList){
        face[0] = match[face[0]];
        face[1] = match[face[1]];
        face[2] = match[face[2]];
    }
    std::cout << "size faceVertexList" << faceVertexList.size() << std::endl;
    std::cout << "size points" << points.size() << std::endl;
    i = 0;
    for(Vector3& p : points){
        std::cout << i++ << " : " << p << std::endl;
    }
    i = 0;
    std::cout << "faces" << std::endl;
    for(std::vector<size_t>& face : faceVertexList){
        std::cout << i++ << " : " << face[0] << ", " << face[1] << ", " << face[2] << std::endl;
    }
    std::cout << "done" << std::endl;
    return makeManifoldSurfaceMeshAndGeometry(faceVertexList, points);
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
mixMeshes(VertexPositionGeometry& geom, VertexPositionGeometry& geomFlat){

    std::cout << "Mixing" << std::endl;
    std::map<size_t, size_t> mapping;
    std::vector<Vector3> pts;
    std::vector<std::vector<size_t>> faces = geom.mesh.getFaceVertexList();

    for(size_t i = 0; i < geom.mesh.nVertices(); ++i)
        pts.push_back(geom.vertexPositions[i]);
    for(Vertex v : geomFlat.mesh.boundaryLoop(0).adjacentVertices()){
        for(Vertex vb : geom.mesh.boundaryLoop(1).adjacentVertices()){
            if(abs(norm(geom.vertexPositions[vb] - geomFlat.vertexPositions[v])) < 10e-8){
                mapping[v.getIndex()] = vb.getIndex();
                break;
            }
        }
    }
    // Add vertices
    for(Vertex v : geomFlat.mesh.vertices()){
        if(mapping.count(v.getIndex()) == 0){
            pts.push_back(geomFlat.vertexPositions[v]);
            std::cout << "added" << geomFlat.vertexPositions[v] << std::endl;
        }
    }

    // Add faces
    size_t n = geom.mesh.nVertices() - mapping.size();
    for(Face f : geomFlat.mesh.faces()){
        size_t a = f.halfedge().vertex().getIndex();
        size_t b = f.halfedge().next().vertex().getIndex();
        size_t c = f.halfedge().next().next().vertex().getIndex();
        if(mapping.count(a) > 0)
            a = mapping[a];
        else
            a += n;
        if(mapping.count(b) > 0)
            b = mapping[b];
        else
            b += n;
        if(mapping.count(c) > 0)
            c = mapping[c];
        else
            c += n;
        faces.push_back({a, b, c});
    }

    for(size_t i = 0; i < pts.size(); ++i)
        std::cout << i << " : " << pts[i] << std::endl;

    for(size_t i = 0; i < faces.size(); ++i)
        std::cout << i << " : " << faces[i][0] << ", " << faces[i][1] << ", " << faces[i][2] << std::endl;

    return makeManifoldSurfaceMeshAndGeometry(faces, pts);
}