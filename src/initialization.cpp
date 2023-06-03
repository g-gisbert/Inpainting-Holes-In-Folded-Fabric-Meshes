#include "initialization.hpp"

void initErrorCorrected(VertexPositionGeometry& geom, BoundaryLoop BL, Plane plane,
                        const std::vector<double>& boundaryEdgeLengths,
                        const std::vector<double>& boundaryAngles) {

    // Get vertex defect angles
    SurfaceMesh& mesh = geom.mesh;
    geom.requireVertexAngleSums();
    std::vector<double> anglesDefect(boundaryAngles);
    int j = 0;
    geom.requireVertexAngleSums();
    for (Vertex v : BL.adjacentVertices()) {
        anglesDefect[j++] = 2*M_PI - geom.vertexAngleSums[v];
    }
    geom.unrequireVertexAngleSums();
    std::vector<double> lengths(boundaryEdgeLengths);
    j = 0;
    for (const Halfedge& he : BL.adjacentHalfedges()) {
        lengths[j++] = norm(geom.vertexPositions[he.vertex()] -
                                geom.vertexPositions[he.next().vertex()]);
    }

    // Move first 2 vertices
    Plane p = Plane(Vector3{0.0f, 1.0f, 0.0f}, plane.shift);
    geom.vertexPositions[BL.halfedge().vertex()] =
            p.projection(geom.vertexPositions[BL.halfedge().vertex()]);
    Vector3 prevPosition = p.projection(geom.vertexPositions[BL.halfedge().vertex()]);
    Vector3 prevDirection = normalize(p.projection(geom.vertexPositions[BL.halfedge().next().vertex()]) -
                                      p.projection(geom.vertexPositions[BL.halfedge().vertex()]));

    double x = cos(M_PI-M_PI/2)*double(prevDirection.x) - sin(M_PI-M_PI/2)*double(prevDirection.z);
    double y = sin(M_PI-M_PI/2)*double(prevDirection.x) + cos(M_PI-M_PI/2)*double(prevDirection.z);
    prevDirection = Vector3{x, 0.0f, y};

    // Project BL vertices
    int i = 0;
    std::vector<Vector3> forwardProject;
    forwardProject.push_back(prevPosition);
    Vector3 error = Vector3::zero();
    for (Vertex v : BL.adjacentVertices()){
        if(i < 1){
            i++;
            continue;
        }

        double theta = -anglesDefect[i-1];
        double l = lengths[i-1];

        double x = cos(M_PI-theta)*double(prevDirection.x) - sin(M_PI-theta)*double(prevDirection.z);
        double y = sin(M_PI-theta)*double(prevDirection.x) + cos(M_PI-theta)*double(prevDirection.z);
        prevDirection = Vector3{x, 0.0f, y}.normalize();

        prevPosition = prevPosition + prevDirection*l;
        forwardProject.push_back(prevDirection*l);
        //geom.vertexPositions[v] = prevPosition;
        error += prevDirection*l;
        i++;
    }

    double theta = -anglesDefect[i-1];
    double l = lengths[0];

    x = cos(M_PI-theta)*double(prevDirection.x) - sin(M_PI-theta)*double(prevDirection.z);
    y = sin(M_PI-theta)*double(prevDirection.x) + cos(M_PI-theta)*double(prevDirection.z);
    prevDirection = Vector3{x, 0.0f, y};
    error += prevDirection*l;

    std::cout << "error : " << error << std::endl;

    i = 0;
    Vector3 lastPos = forwardProject[0];
    for(Vertex v : BL.adjacentVertices()){
        if(i < 1){
            i++;
            continue;}

        geom.vertexPositions[v] = lastPos + forwardProject[i];// - error/BL.degree();
        lastPos = geom.vertexPositions[v];
        i++;
    }

}
