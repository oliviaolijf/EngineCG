#ifndef ENGINE_THREEDEELINES_H
#define ENGINE_THREEDEELINES_H

#include "LSystems.h"
#include "vector/vector3d.h"


class Face{
public:
    std::vector<int> point_indexes;
};

class Figure{
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Color color;
};
typedef std::list<Figure> Figures3D;

#endif //ENGINE_THREEDEELINES_H
