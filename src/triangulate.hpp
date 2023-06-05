#ifndef INC_2DCONTOUR_TRIANGULATE_HPP
#define INC_2DCONTOUR_TRIANGULATE_HPP

#include "pch.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;
typedef CGAL::Spatial_sort_traits_adapter_2<K,Point*> Search_traits;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point3;
typedef CGAL::Surface_mesh<Point3>                           Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor        vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor      halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;


void triangulate(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh, double size);
void triangulateLiepaFair(std::unique_ptr<VertexPositionGeometry>& geom, std::unique_ptr<ManifoldSurfaceMesh>& mesh);


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
toGC(CDT& cdt, std::map<Vertex_handle, std::size_t>& verticesIndex);

template <class InputIterator>
void insert_with_info(CDT& cdt, InputIterator first,InputIterator last)
{
    std::vector<std::ptrdiff_t> indices;
    std::vector<Point> points;
    std::ptrdiff_t index=0;

    for (InputIterator it=first;it!=last;++it){
        points.push_back( *it);
        indices.push_back(index++);
    }

    CGAL::spatial_sort(indices.begin(),indices.end(),Search_traits(&(points[0]),cdt.geom_traits()));

    CDT::Vertex_handle v_hint;
    CDT::Face_handle hint;
    for (typename std::vector<std::ptrdiff_t>::const_iterator
                 it = indices.begin(), end = indices.end();
         it != end; ++it){
        v_hint = cdt.insert(points[*it], hint);
        if (v_hint!=CDT::Vertex_handle()){
            v_hint->info()=*it;
            hint=v_hint->face();
        }
    }
}

#endif //INC_2DCONTOUR_TRIANGULATE_HPP
