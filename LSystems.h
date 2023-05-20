#ifndef ENGINE_LSYSTEMS_H
#define ENGINE_LSYSTEMS_H
#include <iostream>
#include <list>
#include "easy_image.h"
#include "vector/vector3d.h"

class Color {
public:
    double red;
    double blue;
    double green;

    explicit Color() : red(0), green(0), blue(0) {}
    Color(double a, double b, double c){red = a; green = b; blue = c;};
};

class Point2D {
public:
    double x;
    double y;

    explicit Point2D() : x(0), y(0) {}
    Point2D(double a,double b){x = a; y = b;};
};

class Line2D{
public:
    Point2D p1;
    Point2D p2;
    Color color;

    double z1;
    double z2;

    explicit Line2D() : p1(), p2(), color() {}
    Line2D(Point2D a, Point2D b, Color c) {p1 = a; p2 = b; color = c;};


};

using Lines2D = std::list<Line2D>;

#endif //ENGINE_LSYSTEMS_H
