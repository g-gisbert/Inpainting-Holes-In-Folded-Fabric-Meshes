#include "chain.hpp"


double Chain::scale = 0.0;
Vector2 Chain::t = Vector2::zero();

void Chain::saveTargetEdgeLengths() {
    /*edgeLengths.resize(positions.size() - 1);
    for (size_t i = 0; i < positions.size() - 1; ++i) {
        edgeLengths[i] = 1.0*norm(positions[i+1] - positions[i]);
    }*/
    const size_t N = positions.size();
    edgeLengths.resize(N);
    for (size_t i = 0; i < N; ++i) {
        edgeLengths[i] = 1.0*norm(positions[(i+1)%N] - positions[i]);
    }
}

void Chain::saveTargetAngles() {
    const size_t N = positions.size();
    angles.resize(N);
    for (int i = 0; i < int(N); ++i) {
        int id = ((i-1) < 0) ? i-1+N : i-1;
        Vector2 p1 = positions[id];
        Vector2 p2 = positions[i];
        Vector2 p3 = positions[(i+1)%N];

        Vector2 incomingVector = p2 - p1;
        Vector2 outgoingVector = p3 - p2;

        double angle = atan2(-incomingVector.y, -incomingVector.x) - atan2(outgoingVector.y, outgoingVector.x);
        angle = (angle < 0) ? angle+2*M_PI : angle ;
        //angles[i] = angle;
        angles[i] = M_PI;
    }
}

void Chain::setPositions(Vector2 p) {
    positions.push_back(p);
}

void Chain::toSVG() {
    /*Vector2 center;
    for (size_t i = 0; i < positions.size(); ++i) {
        center += ;
    }
    center /= positions.size();*/
    svg::Dimensions dimensions(1000, 1000);
    svg::Document doc("../../svg/my_svg" + std::to_string(frame++) + ".svg", svg::Layout(dimensions, svg::Layout::BottomLeft));

    svg::Polyline polyline(svg::Stroke(2, svg::Color(0,0,0)));

    for (const Vector2& p : positions) {
        svg::Point pt(scale*p.x + t.x, scale*p.y + t.y);
        polyline << pt;
        doc << svg::Circle(pt, 5, svg::Fill(svg::Color(0, 0, 0)));
    }

    doc << polyline;
    bool ok = doc.save();
}

void Chain::toSVG(Eigen::VectorXd& x, Eigen::VectorXd& grad, int k) {
    //std::cout << "ok" << std::endl;
    svg::Dimensions dimensions(1000, 1000);
    svg::Document doc("../../svg/my_svg" + std::to_string(k) + ".svg", svg::Layout(dimensions, svg::Layout::BottomLeft));

    svg::Polyline polyline(svg::Stroke(2, svg::Color(0,0,0)));

    for (size_t i = 0; i < size_t(x.size())/2; ++i) {
        svg::Point pt(scale*x[2*i] + t.x, scale*x[2*i+1] + t.y);
        polyline << pt;
        doc << svg::Circle(pt, 5, svg::Fill(svg::Color(0, 0, 0)));
    }

    // Grad
    constexpr double c = 10000.0;
    for (size_t i = 0; i < size_t(grad.size())/2; ++i) {
        doc << svg::Line(svg::Point(scale*x[2*i] + t.x, scale*x[2*i+1] + t.y),
                         svg::Point(scale*x[2*i]+t.x + c*grad[2*i], scale*x[2*i+1]+t.y + c*grad[2*i+1]),
                         svg::Stroke(2, svg::Color(0,255,0)));
    }

    doc << polyline;
    doc.save();
}

void Chain::readFromMesh(std::string filename) {
    const std::string DATA_PATH = "../../data/";
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH + filename);
    BoundaryLoop BL = mesh->boundaryLoop(0);
    edgeLengths.resize(BL.degree());
    angles.resize(BL.degree());

    int i = 0;
    for (const Halfedge& he : BL.adjacentHalfedges()) {
        //std::cout << "i : " << i << ", " << he.vertex() << std::endl;
        edgeLengths[i++] = 1.0*100*norm(geometry->vertexPositions[he.vertex()] -
                                geometry->vertexPositions[he.next().vertex()]);
    }

    i = 0;
    geometry->requireVertexAngleSums();
    for (Vertex v : BL.adjacentVertices()) {
        //std::cout << "i : " << i << ", " << v << std::endl;
        angles[i++] = 2*M_PI - geometry->vertexAngleSums[v];
        //angles[i++] = geometry->vertexAngleSums[v];
    }
    geometry->unrequireVertexAngleSums();

    initErrorCorrected(*geometry, BL, Plane{Vector3{0, 0, 1}, 1}, edgeLengths, angles);

    i = 0;
    for (Vertex v : BL.adjacentVertices()) {
        setPositions(Vector2{100*geometry->vertexPositions[v].x, -100*geometry->vertexPositions[v].z});
    }
}

void Chain::computeSVGTransform() {

    Vector2 min{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vector2 max{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};

    Vector2 center = Vector2::zero();

    for (const Vector2& p : positions) {
        if (p.x < min.x)
            min.x = p.x;
        if (p.y < min.y)
            min.y = p.y;
        if (p.x > max.x)
            max.x = p.x;
        if (p.y > max.y)
            max.y = p.y;
        center += p;
    }
    center /= double(positions.size());
    double debugscale = std::min(800.0/(max.x - min.x), 800.0/(max.y - min.y));
    Vector2 debugt = Vector2{500, 500} - center;
    scale = std::min(800.0/(max.x - min.x), 800.0/(max.y - min.y));
    t = Vector2{100, 100} - min*scale;

}

void Chain::setPositionsFromMesh(const std::string& filename) {
    const std::string DATA_PATH = "../../data/";
    std::unique_ptr<ManifoldSurfaceMesh> meshTmp;
    std::unique_ptr<VertexPositionGeometry> geomTmp;
    std::tie(meshTmp, geomTmp) = readManifoldSurfaceMesh(DATA_PATH + filename);
    BoundaryLoop BL = meshTmp->boundaryLoop(0);
    int i = 0;
    for (Vertex v : BL.adjacentVertices()) {
        //setPositions(Vector2{100*geometry->vertexPositions[v].x, -100*geometry->vertexPositions[v].z});
        positions[i++] = Vector2{100*geomTmp->vertexPositions[v].x, -100*geomTmp->vertexPositions[v].z};
    }

}
