#ifndef PCH_H
#define PCH_H

// Utilities
#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>
#include <functional>
#include <algorithm>
#include <utility>
#include <tuple>
#include <string>
#include <cstring>
#include <chrono>
#include <fstream>

// STD Containers
#include <queue>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <deque>

// Geometry Central
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/pointcloud/point_cloud.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Polyscope
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

// Eigen
#include <Eigen/Core>

// LibIgl
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/lscm.h>
#include <igl/boundary_loop.h>

#include <igl/arap.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>



// TinyAD
#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <boost/lexical_cast.hpp>

// System
#include <dirent.h>
#include <sys/stat.h>

// Misc
#include <omp.h>



#endif //PCH_H
