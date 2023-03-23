//
// Created by olivi on 3/7/2023.
//

#include "LSystems.h"
#include <cmath>

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