#include "ThreeDeeLines.h"
#include <cmath>
#include <limits>

void toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    r = sqrt((point.x*point.x)+(point.y*point.y)+(point.z*point.z));
    theta = std::atan2(point.y, point.x);
    phi = std::acos(point.z/r);
}

Matrix eyePointTrans(const Vector3D &eyepoint){
    double theta, phi ,r;
    toPolar(eyepoint, theta,phi,r);
    Matrix V;
    V(1,1) = -sin(theta);
    V(1,2) = (-cos(theta)) * cos(phi);
    V(1,3) = cos(theta)*sin(phi);
    V(2,1) = cos(theta);
    V(2,2) = (-sin(theta))* cos(phi);
    V(2,3) = sin(theta)*sin(phi);
    V(3,2) = sin(phi);
    V(3,3) = cos(phi);
    V(4,3) = -r;
    return V;
}

Point2D doProjectionPoint(const Vector3D &point, const double &d){
    auto z = point.z;
    auto xAccent = d*point.x/-z;
    auto yaccent = d*point.y/-z;
    return Point2D(xAccent, yaccent);
};

Lines2D doProjection(const Figures3D &f){
    double d = 1;
    Lines2D lines;
    for(auto& fig :f){
        for (auto face: fig.faces) {
            for (int i = 0; i < face.point_indexes.size() - 1; i++) {
                Line2D line;
                if (i == face.point_indexes.size() - 1) {
                    line.p1 = doProjectionPoint(fig.points[face.point_indexes[i]], d);
                    line.p2 = doProjectionPoint(fig.points[face.point_indexes[0]], d);

                    line.z1 = fig.points[face.point_indexes[i]].z;
                    line.z2 = fig.points[face.point_indexes[0]].z;
                }
                else{
                    line.p1 = doProjectionPoint(fig.points[face.point_indexes[i]], d);
                    line.p2 = doProjectionPoint(fig.points[face.point_indexes[i+1]], d);

                    line.z1 = fig.points[face.point_indexes[i]].z;
                    line.z2 = fig.points[face.point_indexes[i+1]].z;
                }
                line.color = fig.color;
                lines.push_back(line);
            }
        }
    }
    return lines;
}



Matrix scaleFigure(const double scale){
    Matrix S;
    S(1,1) = scale;
    S(2,2) = scale;
    S(3,3) = scale;
    return S;
}

Matrix rotateX(const double angle){
    Matrix Mx;
    Mx(2,2) = cos(angle);
    Mx(2,3) = sin(angle);
    Mx(3,2)= -sin(angle);
    Mx(3,3 ) = cos(angle);
    return Mx;
}

Matrix rotateY(const double angle){
    Matrix My;
    My(1,1) = cos(angle);
    My(1,3) = -sin(angle);
    My(3,1) = sin(angle);
    My(3,3 ) = cos(angle);
    return My;
}

Matrix rotateZ(const double angle){
    Matrix Mz;
    Mz(1,1) = cos(angle);
    Mz(1,2) = sin(angle);
    Mz(2,1) = -sin(angle);
    Mz(2,2)= cos(angle);
    return Mz;
}

Matrix translate(const Vector3D &vector){
    Matrix T;
    T(4,1) = vector.x;
    T(4,2)= vector.y;
    T(4,3) = vector.z;
    return T;
}

void applyTransformation(Figure &f, const Matrix &m){
    for (auto &p: f.points){
        p *= m;
    }
}


