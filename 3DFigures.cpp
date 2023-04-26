#include "ThreeDeeLines.h"

Figure generateCube(){
    Figure cube;
    cube.points = {Vector3D::point(1, -1, -1),
                          Vector3D::point(-1, 1, -1),
                          Vector3D::point(1, 1, 1),
                          Vector3D::point(-1, -1, 1),
                          Vector3D::point(1, 1, -1),
                          Vector3D::point(-1, -1, -1),
                          Vector3D::point(1, -1, 1),
                          Vector3D::point(-1, 1, 1)};
    Face face1({0, 4, 2, 6});
    Face face2({4, 1, 7, 2});
    Face face3({1, 5, 3, 7});
    Face face4({5, 0, 6, 3});
    Face face5({6, 2, 7, 3});
    Face face6({0, 5, 1, 4});
    cube.faces = {face1, face2, face3, face4, face5, face6};
}
return *cubeFigure;
}
}
