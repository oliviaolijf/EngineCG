#ifndef ENGINE_THREEDEELINES_H
#define ENGINE_THREEDEELINES_H

#include "LSystems.h"
#include "vector/vector3d.h"
#include <cmath>


class Face{
public:
    explicit Face() = default;
    Face(std::vector<int> points) {point_indexes = points;};
    std::vector<int> point_indexes;
};

class Figure{
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Color color;
};
typedef std::list<Figure> Figures3D;

void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

Matrix eyePointTrans(const Vector3D &eyepoint);
Matrix scaleFigure(const double scale);
Matrix rotateX(const double angle);
Matrix rotateY(const double angle);
Matrix rotateZ(const double angle);
Matrix translate(const Vector3D &vector);
void applyTransformation(Figure &f, const Matrix &m);

Point2D doProjection(const Vector3D &point, const double &d);
Lines2D doProjection(const Figures3D &f);

#endif //ENGINE_THREEDEELINES_H
