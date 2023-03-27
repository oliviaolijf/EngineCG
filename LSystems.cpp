#include "LSystems.h"
#include <cmath>
#include <tuple>
#include <iostream>
#include <map>

img::Color Color::createColor() {
    img::Color kleur(this->red*255, this->green*255, this->blue*255);
    return kleur;
}

std::list<Point2D> Line2D::get_coordinates() {
    std::list<Point2D> points;
    auto xa = ceil(p1.x); auto ya = ceil(p1.y);
    auto xb = ceil(p2.x); auto yb = ceil(p2.y);
    if (xa > xb){
        auto temp = xa;
        xa = xb;
        xb = temp;
    }
    if (xa == xb){
        for (int i = 0; i < std::max(ya,yb) - std::min(ya, yb); i++){
            auto y = std::min(ya,yb);
            points.push_back(Point2D(xa, y+i));
        }
        return points;
    }
    else if (ya == yb){
        for (int i = 0; i < xb-xa; i++){
            points.push_back(Point2D(xa+i, ya));
        }
        return points;
    }
    double m = (yb-ya)/(xb-xa);
    double xi, yi;
    if (0.0 < m <= 1.0){
        for (int i = 0; i < xb-xa; i++){
            xi = xa + i;
            yi = round(ya+(m*i));
            points.push_back(Point2D(xi,yi));
        }
        return points;
    }
    else if (-1 <= m < 0.0){
        for (int i = 0; i < xb-xa; i++){
            xi = xa + i;
            yi = round(ya+(m*i));
            points.push_back(Point2D(xi,yi));
        }
        return points;
    }
    else if (m > 1){
        for (int i = 0; i < yb-ya; i++){
            yi = ya + i;
            xi = round(xa + (i/m));
            points.push_back(Point2D(xi,yi));
        }
        return points;
    }
    else if (m < -1){
        for (int i = 0; i < ya-yb; i++){
            yi = ya - i;
            xi = round(xa - (i/m));
            points.push_back(Point2D(xi,yi));
        }
        return points;
    }
    return points;
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    auto square_root = (point.x*point.x) + (point.y*point.y) + (point.z * point.z);
    r  = sqrt(square_root);
    theta = std::atan2(point.y, point.x);
    phi = std::acos(r);
}

Matrix eyepointTrans(const Vector3D &eyepoint){
    double r, phi, theta;
    toPolar(eyepoint, theta, phi, r);

    Matrix V;
    V(1,1) = -(sin(theta));
    V(1,2) = -(cos(theta)* cos(phi));
    V(1,3) = cos(theta) * sin(phi);
    V(2,1) = cos(theta);
    V(2,2) = -(sin(theta) * sin(phi));
    V(2,3) = sin(theta)*sin(phi);
    V(3,2) = sin(phi);
    V(3,3) = cos(phi);
    V(4,3) = -r;
    V(4,4) = 1;

    return V;
}


void applyTransformation(Figures3D &figs, const Matrix &m){

}