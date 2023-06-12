#ifndef INC_2DCONTOUR_TRIANGULATE_HPP
#define INC_2DCONTOUR_TRIANGULATE_HPP

#include "pch.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Point Point;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point3;
typedef CGAL::Surface_mesh<Point3>                           Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor        vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor      halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;


void triangulate(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh, double size = 1.41);


#endif //INC_2DCONTOUR_TRIANGULATE_HPP
