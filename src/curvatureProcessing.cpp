#include "curvatureProcessing.hpp"
#include "pch.h"

curvatureProcessing::curvatureProcessing() {}


Vector3 lineSphereIntersection(Vector3 center, float r, Vector3 o, Vector3 u){
    float a = 1.0f;
    float b = 2 * dot(u, o - center);
    float c = pow(norm(o - center), 2) - r*r;
    float delta = b*b - 4*a*c;

    float t = (-b + sqrt(delta)) / (2*a);

    return o + t*u;
}



VertexData<curvatureData> curvatureProcessing::calculateCurvature(VertexPositionGeometry& geometry, float r_coeff){
    /** Calculate the curvature of a mesh for each vertex
     * 
    */
    size_t v_id = 0;
    SurfaceMesh& mesh = geometry.mesh;
    std::unique_ptr<VertexData<double>> maxCurvature = std::make_unique<VertexData<double>>();
    VertexData<curvatureData> data(mesh);
    data[mesh.vertex(v_id)].edgeDebug = EdgeData<float>(mesh);
    data[mesh.vertex(v_id)].edgeDebugDerivative = EdgeData<Vector3>(mesh);
    (*maxCurvature) = VertexData<double>(mesh);
    geometry.requireFaceNormals();

    // Calculate bounding box min and max
    Vector3 xmin = 1000 * Vector3{1.0f, 1.0f, 1.0f};
    Vector3 xmax = 1000 * Vector3{-1.0f, -1.0f, -1.0f};
    for(Vertex v : mesh.vertices()){
        if(geometry.inputVertexPositions[v][0] > xmax[0])
            xmax[0] = geometry.inputVertexPositions[v][0];
        if(geometry.inputVertexPositions[v][0] < xmin[0])
            xmin[0] = geometry.inputVertexPositions[v][0];

        if(geometry.inputVertexPositions[v][1] > xmax[1])
            xmax[1] = geometry.inputVertexPositions[v][1];
        if(geometry.inputVertexPositions[v][1] < xmin[1])
            xmin[1] = geometry.inputVertexPositions[v][1];

        if(geometry.inputVertexPositions[v][2] > xmax[2])
            xmax[2] = geometry.inputVertexPositions[v][2];
        if(geometry.inputVertexPositions[v][2] < xmin[2])
            xmin[2] = geometry.inputVertexPositions[v][2];
        data[v].boundary = 0;
    }
    float r = r_coeff * norm(xmax - xmin);
    //std::cout << "xmax : " << xmax << ", xmin : " << xmin << " , r : " << r << std::endl;

    int i = 0;
    for(Vertex v : mesh.vertices()){
        Vector3 v_pos = geometry.vertexPositions[v];
        //std::cout << std::endl << i++ << "/" << mesh.nVertices() << std::endl;

        // Retrieve edges in the sphere of radius r
        std::vector<std::pair<Edge, float>> edgesInSphere;
        std::map<int, bool> verticesInSphere;
        for(Edge e : mesh.edges()){
            Vertex v1 = e.firstVertex();
            Vertex v2 = e.secondVertex();

            if(norm(geometry.inputVertexPositions[v1] - v_pos) < r){
                verticesInSphere.insert(std::pair<int, bool>(v1.getIndex(), true));
                if(v.getIndex() == v_id) data[v1].boundary = 1;
                if(norm(geometry.inputVertexPositions[v2] - v_pos) < r){
                    verticesInSphere.insert(std::pair<int, bool>(v2.getIndex(), true));
                    edgesInSphere.push_back(std::pair<Edge, float>{e, 1.0f});
                }
                else{
                    edgesInSphere.push_back(std::pair<Edge, float>{e, 0.5f});

                }

            }
            else if(norm(geometry.inputVertexPositions[v2] - v_pos) < r){
                verticesInSphere.insert(std::pair<int, bool>(v2.getIndex(), true));
                edgesInSphere.push_back(std::pair<Edge, float>{e, 0.5f});
                if(v.getIndex() == v_id) data[v2].boundary = 1;
            }
        }
        float B = 0.0f;
        for(const auto &pair : verticesInSphere){
            int v1_id = pair.first;
            B += geometry.vertexDualArea(mesh.vertex(v1_id));
        }

        // Computes d|B|/dvi
        Vector3 dB = Vector3{0.0f, 0.0f, 0.0f};
        std::map<int, bool> facesInSphere;
        for(std::pair<Edge, float>& p : edgesInSphere){
            Edge e = p.first;
            int f1_id = e.halfedge().face().getIndex();
            int f2_id = e.halfedge().twin().face().getIndex();
            facesInSphere.insert(std::pair<int, bool>(f1_id, true));
            if(e.halfedge().twin().isInterior())
                facesInSphere.insert(std::pair<int, bool>(f2_id, true));
        }
        for ( const auto &pair : facesInSphere ) {
            int f_id = pair.first;
            Face f = mesh.face(f_id);
            Vertex v1 = f.halfedge().vertex();
            Vertex v2 = f.halfedge().next().vertex();
            Vertex v3 = f.halfedge().next().next().vertex();

            // Bf
            float Bf = 0.0f;
            if(norm(geometry.inputVertexPositions[v1] - v_pos) < r)
                Bf += 0.33333333f;
            if(norm(geometry.inputVertexPositions[v2] - v_pos) < r)
                Bf += 0.33333333f;
            if(norm(geometry.inputVertexPositions[v3] - v_pos) < r)
                Bf += 0.33333333f;

            // dAf
            Vector3 n = geometry.faceNormal(f);
            Halfedge heCurrent = f.halfedge();

            bool isInside = false;
            for(int k = 0; k < 3; ++k){
                if(v.getIndex() != heCurrent.vertex().getIndex())
                    heCurrent = heCurrent.next();
                else
                    isInside = true;
            }
            if(isInside){
                Vector3 e23 = geometry.vertexPositions[heCurrent.next().next().vertex().getIndex()] -
                              geometry.vertexPositions[heCurrent.next().vertex().getIndex()];

                dB += Bf * 0.5*cross(n, e23);
            }

        }

        Eigen::Matrix3f K;
        K << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
        Eigen::Matrix3f dKx;
        dKx << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
        Eigen::Matrix3f dKy;
        dKy << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
        Eigen::Matrix3f dKz;
        dKz << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;

        Eigen::Matrix3f dKx_hat;
        dKx_hat << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
        Eigen::Matrix3f dKy_hat;
        dKy_hat << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
        Eigen::Matrix3f dKz_hat;
        dKz_hat << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;


        // for e in B
        for(std::pair<Edge, float>& p : edgesInSphere){
            Edge e = p.first;
            float eInterB = p.second;

            // Check if v is v1 or v2
            Vertex v1 = e.firstVertex();
            Vertex v2 = e.secondVertex();
            Vector3 v1_pos = geometry.vertexPositions[v1];
            Vector3 v2_pos = geometry.vertexPositions[v2];
            int whichIsV = -1;
            if(v.getIndex() == v1.getIndex())
                whichIsV = 1;
            if(v.getIndex() == v2.getIndex())
                whichIsV = 2;


            // Edge direction
            Vector3 eLength = v1_pos - v2_pos;
            float eNorm = norm(eLength);
            Vector3 eDirection = eLength / eNorm;
            eInterB = norm(eInterB*eLength);

            // Matrix construction
            K(0,0) += eDirection[0]*eDirection[0]*geometry.edgeDihedralAngle(e)*eInterB;K(1,0) += eDirection[1]*eDirection[0]*geometry.edgeDihedralAngle(e)*eInterB;
            K(2,0) += eDirection[2]*eDirection[0]*geometry.edgeDihedralAngle(e)*eInterB;
            K(0,1) += eDirection[0]*eDirection[1]*geometry.edgeDihedralAngle(e)*eInterB;K(1,1) += eDirection[1]*eDirection[1]*geometry.edgeDihedralAngle(e)*eInterB;
            K(2,1) += eDirection[2]*eDirection[1]*geometry.edgeDihedralAngle(e)*eInterB;
            K(0,2) += eDirection[0]*eDirection[2]*geometry.edgeDihedralAngle(e)*eInterB;K(1,2) += eDirection[1]*eDirection[2]*geometry.edgeDihedralAngle(e)*eInterB;
            K(2,2) += eDirection[2]*eDirection[2]*geometry.edgeDihedralAngle(e)*eInterB;


            Eigen::Matrix3f eeT_bar;
            eeT_bar << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
            eeT_bar(0,0) = eDirection[0]*eDirection[0];eeT_bar(1,0) = eDirection[1]*eDirection[0];
            eeT_bar(2,0) = eDirection[2]*eDirection[0];
            eeT_bar(0,1) = eDirection[0]*eDirection[1];eeT_bar(1,1) = eDirection[1]*eDirection[1];
            eeT_bar(2,1) = eDirection[2]*eDirection[1];
            eeT_bar(0,2) = eDirection[0]*eDirection[2];eeT_bar(1,2) = eDirection[1]*eDirection[2];
            eeT_bar(2,2) = eDirection[2]*eDirection[2];

            // deInterB
            Vector3 deInterB = Vector3{0.0f, 0.0f, 0.0f};
            float Be12 = 0.5f;
            if(norm(geometry.vertexPositions[v1] - v_pos) < r && norm(geometry.vertexPositions[v2] - v_pos) < r)
                Be12 = 1.0f;

            if(whichIsV == 1){
                deInterB = Be12 * eDirection;
            }
            else if(whichIsV == 2){
                deInterB = -Be12 * eDirection;
            }


            // dBeta
            Vector3 dBeta = Vector3{0.0f, 0.0f, 0.0f};
            Face f0 = e.halfedge().face();
            Face f1 = e.halfedge().twin().face();
            if(e.halfedge().twin().isInterior()){
                Vector3 n0 = geometry.faceNormal(f0);
                Vector3 n1 = geometry.faceNormal(f1);
                if(v.getIndex() == e.halfedge().vertex().getIndex()){
                    dBeta = 1/eNorm * (2*geometry.halfedgeCotanWeight(e.halfedge().next().next())*n0 +
                                       2*geometry.halfedgeCotanWeight(e.halfedge().twin().next())*n1);
                }
                if(v.getIndex() == e.halfedge().twin().vertex().getIndex()){
                    dBeta = 1/eNorm * (2*geometry.halfedgeCotanWeight(e.halfedge().next())*n0 +
                                       2*geometry.halfedgeCotanWeight(e.halfedge().twin().next().next())*n1);
                }
                if(v.getIndex() == e.halfedge().next().next().vertex().getIndex()){
                    //dBeta = -eNorm*n0/(2*geometry.faceArea(f0));
                }
                if(v.getIndex() == e.halfedge().twin().next().next().vertex().getIndex()){
                    //dBeta = -eNorm*n1/(2*geometry.faceArea(f1));
                }
            }



            // deeT_bar
            Eigen::Matrix3f deeT_bar_x; deeT_bar_x << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
            Eigen::Matrix3f deeT_bar_y; deeT_bar_y << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
            Eigen::Matrix3f deeT_bar_z; deeT_bar_z << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;

            Eigen::Matrix3f deeT_x; deeT_x << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
            Eigen::Matrix3f deeT_y; deeT_y << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
            Eigen::Matrix3f deeT_z; deeT_z << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;

            Vector3 de12_norm = Vector3{0.0f, 0.0f, 0.0f};

            deeT_x(0,0) = 2*eLength[0]; deeT_x(0,1) = eLength[1];
            deeT_x(0,2) = eLength[2];
            deeT_x(1,0) = eLength[1]; deeT_x(1,1) = 0.0f;
            deeT_x(1,2) = 0.0f;
            deeT_x(2,0) = eLength[2]; deeT_x(2,1) = 0.0f;
            deeT_x(2,2) = 0.0f;

            deeT_y(0,0) = 0.0f; deeT_y(0,1) = eLength[0];
            deeT_y(0,2) = 0.0f;
            deeT_y(1,0) = eLength[0]; deeT_y(1,1) = 2*eLength[1];
            deeT_y(1,2) = eLength[2];
            deeT_y(2,0) = 0.0f; deeT_y(2,1) = eLength[2];
            deeT_y(2,2) = 0.0f;

            deeT_z(0,0) = 0.0f; deeT_z(0,1) = 0.0f;
            deeT_z(0,2) = eLength[0];
            deeT_z(1,0) = 0.0f; deeT_z(1,1) = 0.0f;
            deeT_z(1,2) = eLength[1];
            deeT_z(2,0) = eLength[0]; deeT_z(2,1) = eLength[1];
            deeT_z(2,2) = 2*eLength[2];

            if(whichIsV == 1){
                de12_norm = 2 * eLength;
            }
            else if (whichIsV == 2){
                de12_norm = -2 * eLength;
                deeT_x = -deeT_x;
                deeT_y = -deeT_y;
                deeT_z = -deeT_z;
            }
            else{
                deeT_x << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
                deeT_y << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
                deeT_z << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
            }


            deeT_bar_x = deeT_x / (eNorm*eNorm) - eeT_bar / (eNorm*eNorm) * de12_norm.x;
            deeT_bar_y = deeT_y / (eNorm*eNorm) - eeT_bar / (eNorm*eNorm) * de12_norm.y;
            deeT_bar_z = deeT_z / (eNorm*eNorm) - eeT_bar / (eNorm*eNorm) * de12_norm.z;

            if(v.getIndex() == v_id){
                data[v].edgeDebugDerivative[e] = dBeta;
                data[v].edgeDebug[e] = geometry.edgeDihedralAngle(e);
            }

            dKx_hat += dBeta.x * eInterB * eeT_bar +
                       geometry.edgeDihedralAngle(e) * deInterB.x * eeT_bar +
                       geometry.edgeDihedralAngle(e) * eInterB * deeT_bar_x;
            dKy_hat += dBeta.y * eInterB * eeT_bar +
                       geometry.edgeDihedralAngle(e) * deInterB.y * eeT_bar +
                       geometry.edgeDihedralAngle(e) * eInterB * deeT_bar_y;
            dKz_hat += dBeta.z * eInterB * eeT_bar +
                       geometry.edgeDihedralAngle(e) * deInterB.z * eeT_bar +
                       geometry.edgeDihedralAngle(e) * eInterB * deeT_bar_z;

            /*if(i == 4){
            std::cout << "########## i : " << i << std::endl;
            std::cout << "e : " << e << std::endl;
            
            std::cout << "deeT_x : " << deeT_x << std::endl;
            std::cout << "deeT_y : " << deeT_y << std::endl;
            std::cout << "deeT_z : " << deeT_z << std::endl;
            std::cout << "eNorm : " << eNorm << std::endl;
            std::cout << "eLength : " << eLength << std::endl;
            std::cout << "de12_norm : " << de12_norm << std::endl;
            std::cout << "dBeta : " << dBeta << std::endl;
            std::cout << "beta : " << geometry.edgeDihedralAngle(e) << std::endl;
            std::cout << "eInterB : " << eInterB << std::endl;
            std::cout << "deInterB : " << deInterB << std::endl;
            std::cout << "eeT_bar : " << eeT_bar << std::endl;
            std::cout << "deeT_bar_x : " << deeT_bar_x << std::endl;
            std::cout << "deeT_bar_y : " << deeT_bar_y << std::endl;
            std::cout << "deeT_bar_z : " << deeT_bar_z << std::endl;
            std::cout << "term_x : " << dBeta.x * eInterB * eeT_bar << std::endl;
            std::cout << "term2_x : " << geometry.edgeDihedralAngle(e) * deInterB.x * eeT_bar << std::endl;
            std::cout << "term3_x : " << geometry.edgeDihedralAngle(e) * eInterB * deeT_bar_x << std::endl;
            std::cout << "term_y : " << dBeta.y * eInterB * eeT_bar << std::endl;
            std::cout << "term2_y : " << geometry.edgeDihedralAngle(e) * deInterB.y * eeT_bar << std::endl;
            std::cout << "term3_y : " << geometry.edgeDihedralAngle(e) * eInterB * deeT_bar_y << std::endl;
            std::cout << "term_z : " << dBeta.z * eInterB * eeT_bar << std::endl;
            std::cout << "term2_z : " << geometry.edgeDihedralAngle(e) * deInterB.z * eeT_bar << std::endl;
            std::cout << "term3_z : " << geometry.edgeDihedralAngle(e) * eInterB * deeT_bar_z << std::endl;
            }*/

        }
        Eigen::Matrix3f K_hat = K;
        K /= B;
        //B = geometry.vertexDualArea(v);
        dKx = 1/B * dKx_hat - (K_hat)/(B*B) * dB.x;
        dKy = 1/B * dKy_hat - (K_hat)/(B*B) * dB.y;
        dKz = 1/B * dKz_hat - (K_hat)/(B*B) * dB.z;


        /*if(i == 9) std::cout << "dKx_hat : " << dKx_hat << std::endl;
        if(i == 9) std::cout << "dB : " << dB << std::endl;
        if(i == 9) std::cout << "dKx : " << dKx << std::endl;*/

        /*if(i == 0){
            for(Edge e : mesh.edges())
                data[0].edgeDebug[e] = geometry.edgeDihedralAngle(e);
        }*/



        /*if(i == 4){
            std::cout << std::endl << "######" << std::endl;
            std::cout << "K : " << K << std::endl;
            std::cout << "edgesInSphere : " << edgesInSphere.size() << std::endl;
            std::cout << "dKx : " << dKx << std::endl;
            std::cout << "dKx_hat : " << dKx_hat << std::endl;
            std::cout << "dKy : " << dKy << std::endl;
            std::cout << "dKy_hat : " << dKy_hat << std::endl;
            std::cout << "dKz : " << dKz << std::endl;
            std::cout << "dKz_hat : " << dKz_hat << std::endl;
            std::cout << "B : " << B << std::endl;
            std::cout << "dB : " << dB << std::endl;
        }*/
        i++;

        // Computes k1 k2  
        Eigen::EigenSolver<Eigen::Matrix3f> es(K); // eigen values and vectors solver
        Eigen::Vector3f eigenValues = es.eigenvalues().real();
        /*std::cout << "EigenVector 1 : " << es.eigenvectors().col(0).real() << std::endl;
        std::cout << "EigenVector 2 : " << es.eigenvectors().col(1).real() << std::endl;
        std::cout << "EigenVector 3 : " << es.eigenvectors().col(2).real() << std::endl;
        if(i > 2000){
            std::cout << "Eigenvalues : " <<  es.eigenvalues() << std::endl;
        }*/

        float min = std::numeric_limits<float>::max(); int id = 0;
        for(int j = 0; j < 3; ++j){
            if(abs(eigenValues.real()[j]) < min){
                min = abs(eigenValues.real()[j]);
                id = j;
            }
        }
        float kn = eigenValues[id];
        data[v].n = Vector3{es.eigenvectors().col(id).real()[0],
                            es.eigenvectors().col(id).real()[1],
                            es.eigenvectors().col(id).real()[2]};
        float k1k2[2] = {0.0f, 0.0f};
        Vector3 u1u2[2] = {Vector3{0.0f, 0.0f, 0.0f}, Vector3{0.0f, 0.0f, 0.0f}};
        int cpt = 0;
        for(int j = 0; j < 3; ++j){
            if(j != id){
                k1k2[cpt] = eigenValues.real()[j];
                u1u2[cpt++] = Vector3{es.eigenvectors().col(j).real()[0],
                                      es.eigenvectors().col(j).real()[1],
                                      es.eigenvectors().col(j).real()[2]};
            }
        }

        if(abs(k1k2[0]) > abs(k1k2[1])){
            data[v].k1 = k1k2[0];
            data[v].u1 = u1u2[0];
            data[v].k2 = k1k2[1];
            data[v].u2 = u1u2[1];
        }
        else{
            data[v].k1 = k1k2[1];
            data[v].u1 = u1u2[1];
            data[v].k2 = k1k2[0];
            data[v].u2 = u1u2[0];
        }

        // Computes derivatives
        float epsilon = 10e-2;
        // Case 3 - umbilic point
        if(abs(data[v].k1 - data[v].k2) < epsilon){
            //std::cout << "umbilic point" << std::endl;
            data[v].dk1x = 0.5f * (dKx(0,0) + dKx(1,1) + dKx(2,2));
            data[v].dk1y = 0.5f * (dKy(0,0) + dKy(1,1) + dKy(2,2));
            data[v].dk1z = 0.5f * (dKz(0,0) + dKz(1,1) + dKz(2,2));
            data[v].dk2x = data[v].dk1x;
            data[v].dk2y = data[v].dk1y;
            data[v].dk2z = data[v].dk1z;

            data[v].pointType = 3;
        }


            // Case 2 - cylindric point
        else if (abs(data[v].k1 - kn) < epsilon && abs(data[v].k1 - data[v].k2) > 3*epsilon){
            //std::cout << "cylindric point k1" << std::endl;
            Eigen::Vector3f u2(data[v].u2.x, data[v].u2.y, data[v].u2.z);

            data[v].dk2x = u2.transpose() * dKx * u2;
            data[v].dk2y = u2.transpose() * dKy * u2;
            data[v].dk2z = u2.transpose() * dKz * u2;

            data[v].dk1x = dKx(0,0) + dKx(1,1) + dKx(2,2) - data[v].dk2x;
            data[v].dk1y = dKy(0,0) + dKy(1,1) + dKy(2,2) - data[v].dk2y;
            data[v].dk1z = dKz(0,0) + dKz(1,1) + dKz(2,2) - data[v].dk2z;

            data[v].pointType = 2;
        }
        else if(abs(data[v].k2 - kn) < epsilon && abs(data[v].k1 - data[v].k2) > 3*epsilon){
            //std::cout << "cylindric point k2" << std::endl;
            Eigen::Vector3f u1(data[v].u1.x, data[v].u1.y, data[v].u1.z);

            data[v].dk1x = u1.transpose() * dKx * u1;
            data[v].dk1y = u1.transpose() * dKy * u1;
            data[v].dk1z = u1.transpose() * dKz * u1;

            data[v].dk2x = dKx(0,0) + dKx(1,1) + dKx(2,2) - data[v].dk1x;
            data[v].dk2y = dKy(0,0) + dKy(1,1) + dKy(2,2) - data[v].dk1y;
            data[v].dk2z = dKz(0,0) + dKz(1,1) + dKz(2,2) - data[v].dk1z;

            data[v].pointType = 2;
        }
            // Case 1 - general case
        else {
            //std::cout << "general case" << std::endl;
            Eigen::Vector3f u1(data[v].u1.x, data[v].u1.y, data[v].u1.z);
            Eigen::Vector3f u2(data[v].u2.x, data[v].u2.y, data[v].u2.z);

            data[v].dk1x = u1.transpose() * dKx * u1;
            data[v].dk1y = u1.transpose() * dKy * u1;
            data[v].dk1z = u1.transpose() * dKz * u1;

            data[v].dk2x = u2.transpose() * dKx * u2;
            data[v].dk2y = u2.transpose() * dKy * u2;
            data[v].dk2z = u2.transpose() * dKz * u2;

            data[v].pointType = 1;
        }
        data[v].dHx = 0.5f * (dKx(0,0) + dKx(1,1) + dKx(2,2));
        data[v].dHy = 0.5f * (dKy(0,0) + dKy(1,1) + dKy(2,2));
        data[v].dHz = 0.5f * (dKz(0,0) + dKz(1,1) + dKz(2,2));
        data[v].H = 0.5f*(data[v].k1 + data[v].k2);
    }

    return data;
}






/*
Previous versions :

// Signed angle
Face faces[2]; int cpt = 0;
for(Face f : e.adjacentFaces()){
    faces[cpt++] = f;
}
float angle = std::acos(dot( geometry.faceNormals[faces[0]], geometry.faceNormals[faces[1]] ));
Vector3 rotatedNormal = cross(geometry.faceNormals[faces[0]], eDirection);
if(dot(rotatedNormal, geometry.faceNormals[faces[1]]) < 0)   
    angle = -angle;
//std::cout << "norm " << angle << std::endl;
//std::cout << "eInterB " << eInterB << std::endl;


// Retrieve edges in the sphere of radius r
std::vector<std::pair<Edge, float>> edgesInSphere;
for(Edge e : mesh.edges()){
    Vertex v1 = e.firstVertex();
    Vertex v2 = e.secondVertex();

    if(norm(geometry.inputVertexPositions[v1] - v_pos) < r){

        Vector3 u = geometry.vertexPositions[v2] - geometry.vertexPositions[v1];

        if(norm(geometry.inputVertexPositions[v2] - v_pos) < r){
            float length = norm(u);
            edgesInSphere.push_back(std::pair<Edge, float>{e, length});
        }
        else{
            Vector3 P = lineSphereIntersection(v_pos, r, geometry.inputVertexPositions[v1], u/norm(u));
            float length = norm(P - geometry.inputVertexPositions[v1]);
            edgesInSphere.push_back(std::pair<Edge, float>{e, length});
        }
    }
    else if(norm(geometry.inputVertexPositions[v2] - v_pos) < r){
        Vector3 u = geometry.vertexPositions[v1] - geometry.vertexPositions[v2];
        Vector3 P = lineSphereIntersection(v_pos, r, geometry.inputVertexPositions[v2], u/norm(u));
        float length = norm(P - geometry.inputVertexPositions[v2]);
        edgesInSphere.push_back(std::pair<Edge, float>{e, length});
    }
}


*/