#include "ThreeDeeLines.h"
#include <vector>
#include <cmath>

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
    V(4,4) = 1;
    V(1,4) = 0;
    V(2, 4) = 0;
    V(3,4) = 0;
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
        auto faces = fig.faces;
        for (auto face: faces){
            int point1 = face.point_indexes[0];
            int point2 = face.point_indexes[1];

            Vector3D vec1 = Vector3D::point(fig.points[point1]);
            Vector3D vec2 = Vector3D::point(fig.points[point2]);

            auto p1 = doProjectionPoint(vec1, d);
            auto p2 = doProjectionPoint(vec2, d);

            Line2D l(p1, p2, fig.color);
            lines.push_back(l);
        }
    }
    return lines;
}
