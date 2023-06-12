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

    PMP::extract_boundary_cycles(meshCGAL, std::back_inserter(border_cycles));
    halfedge_descriptor h = border_cycles[0];
    std::vector<face_descriptor>  patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    PMP::triangulate_and_refine_hole(meshCGAL,
        h, std::back_inserter(patch_facets), std::back_inserter(patch_vertices),
        CGAL::parameters::density_control_factor(size));
    CGAL::IO::write_polygon_mesh("tmp2.obj", meshCGAL, CGAL::parameters::stream_precision(17));
    std::tie(mesh, geom) = readManifoldSurfaceMesh("tmp2.obj");
}
