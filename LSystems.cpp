//
// Created by olivi on 3/7/2023.
//

#include "LSystems.h"
#include <cmath>

std::list<Point2D> Line2D::get_coordinates() {
    double x1 = p1.x, y1 = p1.y;
    double x2 = p2.x, y2 = p2.y;
    double x,y;
    std::list<Point2D> points;
    if (x1 == x2){//verticale lijn
        Point2D p(0,0);
        p.x = x1;
        if (y1 < y2){
            for (int i = 0; i <= (y2-y1); i++){
                y = y1+i;
                p.y = y;
                points.push_back(p);
            }
        }
        else if (y2 < y1){
            for (int i = 0; i <= (y1-y2); i++){
                y = y2+i;
                p.y = y;
                points.push_back(p);
            }
        }
        return points;
    }
    else if (y1 == y2){
        Point2D p(0,0);
        p.y = y1;
        if (x1 < x2){
            for (int i = 0; i <= (x2-x1); i++){
                x = x1+i;
                p.x = x;
                points.push_back(p);
            }
        }
        else if (x2 < x1){
            for (int i = 0; i <= (x1-x2); i++){
                x = x2+i;
                p.x = x;
                points.push_back(p);
            }
        }
        return points;
    }
    double m = (y2-y1)/(x2-x1);
    if (m <= 1 && 0 < m){
        double xi, yi;
        for (int i = 0; i < x2-x1; i++){
            xi = x1+i;
            yi = round(y1 + m*i);
            points.push_back(Point2D(xi,yi));
        }
        return points;
    }
    else if (-1 <= m && m < 0){
        double xi, yi;
        for (int i = 0; i < x2-x1; i++){
            xi = x1+i;
            yi = round(y1 + m*i);
            points.push_back(Point2D(xi,yi));
        }
        return points;
    }
    else if (m > 1) {
        double xi, yi;
        for (int i = 0; i < y2-y1; i++){
            yi = y1+i;
            xi = round(x1+(i/m));
            points.push_back(Point2D(xi,yi));
        }
    }
    else if (m < -1){
        double xi, yi;
        for (int i = 0; i < y1-y2; i++){
            yi = y1+i;
            xi = round(x1+(i/m));
            points.push_back(Point2D(xi,yi));
        }
    }
    return points;
}

std::list<Point2D> Line2D::get_begin_and_end() {
    std::list<Point2D> points;
    points.push_back(p1);
    points.push_back(p2);
    return points;
}
