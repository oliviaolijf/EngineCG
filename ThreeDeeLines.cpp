#include "ThreeDeeLines.h"

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
    if (point.z ==0){
        z = 1;
    }
    auto xAccent = d*point.x/-z;
    auto yaccent = d*point.y/-z;
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

Matrix scaleFigure(const double scale){
    Matrix S;
    S(1,1) = scale;
    S(2,2) = scale;
    S(3,3) = scale;
    return S;
}

Matrix rotateX(const double angle){
    Matrix Mx;
    Mx(1,1) = 1;
    Mx(2,2) = std::cos(angle);
    Mx(2,3) = std::sin(angle);
    Mx(3,2)= -(std::sin(angle));
    Mx(3,3 ) = std::cos(angle);
    return Mx;
}

Matrix rotateZ(const double angle){
    Matrix My;
    My(1,1) = cos(angle);
    My(1,3) = -sin(angle);
    My(2,2)= 1;
    My(3,1) = sin(angle);
    My(3,3 ) = cos(angle);

    return My;
}

Matrix RotateZ(const double angle){
    Matrix Mz;
    Mz(1,1) = cos(angle);
    Mz(1,2) = sin(angle);
    Mz(2,1) = -sin(angle);
    Mz(2,2)= cos(angle);
    Mz(3,3) = 1;

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

}
