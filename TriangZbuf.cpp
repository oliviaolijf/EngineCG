#include "ThreeDeeLines.h"
#include "Zbuffer.h"
#include "LSystems.h"
#include "3DFigures.cpp"

std::vector<Face> triangulate(const Face& face){
    auto n = face.point_indexes.size();
    std::vector<Face> triangles;

    for (int i = 2; i < n-2; i++){
        Face f = {{face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]}};
        triangles.push_back(f);
    }
    return triangles;
}