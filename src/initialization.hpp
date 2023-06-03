#ifndef INC_2DCONTOUR_INITIALIZATION_HPP
#define INC_2DCONTOUR_INITIALIZATION_HPP

#include "pch.h"
#include "plane.hpp"

void initErrorCorrected(VertexPositionGeometry& geom, BoundaryLoop BL, Plane plane,
                        const std::vector<double>& boundaryEdgeLengths,
                        const std::vector<double>& boundaryAngles);

#endif //INC_2DCONTOUR_INITIALIZATION_HPP
