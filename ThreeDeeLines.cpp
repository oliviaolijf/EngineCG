#include "ThreeDeeLines.h"

void toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    r = sqrt((point.x*point.x)+(point.y*point.y)+(point.z*point.z));
    theta = std::atan2(point.y, point.x);
    phi = std::acos(r);
}

Matrix eyePointTrans(const Vector3D &eyepoint){
    double theta, phi ,r;
    toPolar(eyepoint, theta,phi,r);
    Matrix V;
    V(1,1) = -(sin(theta));
    V(2,1) = cos(theta);
    V(1,2) = -(cos(theta)) * cos(phi);
    V(2,2) = -sin(theta)* cos(phi);
    V(3,2) = sin(phi);
    V(1,3) = cos(theta)*sin(phi);
    V(2,3) = sin(theta)*sin(phi);
    V(3,3) = cos(phi);
    V(4,3) = -r;
    return V;
}

Point2D doProjectionPoint(const Vector3D &point, const double &d){
    auto xAccent = (d*point.x)/-(point.z);
    auto yaccent = (d*point.y)/-(point.z);
    return Point2D(xAccent, yaccent);
};

Lines2D doProjection(const Figures3D &f){
    double d = 1;
    Lines2D lines;
    for(auto& fig :f){
        for (auto face: fig.faces){
            Line2D line;
            line.p1 = doProjectionPoint(fig.points[face.point_indexes[0]], d);
            line.p2 = doProjectionPoint(fig.points[face.point_indexes[1]], d);
            lines.push_back(line);
        }
    }
    return lines;
}
