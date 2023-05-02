#include "ThreeDeeLines.h"

Figure generateCube(Color c){
    Figure cube;
    cube.points = {Vector3D::point(1,-1,-1),
                   Vector3D::point(-1,1,-1),
                   Vector3D::point(1,1,1),
                   Vector3D::point(-1,-1,1),
                   Vector3D::point(1,1,-1),
                   Vector3D::point(-1,-1,-1),
                   Vector3D::point(1,-1,1),
                   Vector3D::point(-1,1,1)};
    Face f1({0,4,2,6});
    Face f2({4,1,7,2});
    Face f3({1,5,3,7});
    Face f4({5,2,6,3});
    Face f5({6,2,7,3});
    Face f6({0,4,1,5});
    cube.faces = {f1,f2,f3,f4,f5,f6};
    cube.color = c;
    return cube;
}
