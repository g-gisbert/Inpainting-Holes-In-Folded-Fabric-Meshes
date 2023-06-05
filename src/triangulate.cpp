#include "triangulate.hpp"



void triangulate(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh, double size) {

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
        CGAL::parameters::density_control_factor(1.7));
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
