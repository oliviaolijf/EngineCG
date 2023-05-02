#include "ThreeDeeLines.h"
#include <cmath>
#define M_PI 3.14159265358979


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
    Face f4({5,0,6,3});
    Face f5({6,2,7,3});
    Face f6({0,5,1,4});
    cube.faces = {f1,f2,f3,f4,f5,f6};
    cube.color = c;
    return cube;
}

Figure generateTetrahedron(Color c){
    Figure tetra;
    tetra.points = {Vector3D::point(1,-1,-1),
                    Vector3D::point(-1,1,-1),
                    Vector3D::point(1,1,1),
                    Vector3D::point(-1,-1,1)
                    };
    Face f1({0,1,2});
    Face f2({1,3,2});
    Face f3({0,3,1});
    Face f4({0,2,3});

    tetra.faces = {f1,f2,f3,f4};
    tetra.color = c;
    return tetra;
}

Figure generateOctahedron(Color c){
    Figure octa;
    octa.points = {Vector3D::point(1,0,0),
                   Vector3D::point(0,1,0),
                   Vector3D::point(-1,0,0),
                   Vector3D::point(0,-1,0),
                   Vector3D::point(0,0,-1),
                   Vector3D::point(0,0,1)
                   };
    Face f1({0, 1, 5});
    Face f2({1,2,5});
    Face f3({2,3,5});
    Face f4({3,0,5});
    Face f5({1,0,4});
    Face f6({2,1,4});
    Face f7({3,2,4});
    Face f8({0,3,4});

    octa.faces = {f1,f2,f3,f4,f5,f6,f7,f8};
    octa.color = c;
    return octa;
}

Figure generateIcosahedron(Color c){
    Figure ico;
    ico.points.push_back(Vector3D::point(0,0, sqrt(5)/2));

    for (int i = 2; i <= 6; i++){
        double x = std::cos((i-2)*2*M_PI/5);
        double y = std::sin((i-2)*2*M_PI/5);
        double z = 0.5;

        ico.points.push_back(Vector3D::point(x,y,z));
    }

    for (int i = 7; i <= 11; i++){
        double x = std::cos(M_PI/5 + (i-7)*2*M_PI/5);
        double y = std::sin(M_PI/5 + (i-7)*2*M_PI/5);
        double z = -0.5;

        ico.points.push_back(Vector3D::point(x,y,z));
    }

    ico.points.push_back(Vector3D::point(0,0,-sqrt(5)/2));

    Face f1({0,1,2});
    Face f2({0,2,3});
    Face f3({0,3,4});
    Face f4({0,4,5});
    Face f5({0,5,1});

    Face f6({1,6,2});
    Face f7({2,6,7});
    Face f8({2,7,3});
    Face f9({3,7,8});
    Face f10({3,8,4});

    Face f11({4,8,9});
    Face f12({4,9,5});
    Face f13({5,9,10});
    Face f14({5,10,1});
    Face f15({1,10,6});

    Face f16({11,7,6});
    Face f17({11,8,7});
    Face f18({11,9,8});
    Face f19({11,10,9});
    Face f20({11,6,10});

    ico.faces = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                 f11,f12,f13,f14,f15,f16,f17,f18,f19,f20};

    ico.color = c;
    return ico;
}

Figure generateDodecahedron(Color c){
    Figure Dode;
    auto ico = generateIcosahedron(c);
    double x1  = (ico.points[0].x + ico.points[1].x + ico.points[2].x)/3;
    double y1  = (ico.points[0].y + ico.points[1].y + ico.points[2].y)/3;
    double z1  = (ico.points[0].z + ico.points[1].z + ico.points[2].z)/3;
    Dode.points.push_back(Vector3D::point(x1,y1,z1));

    double x2  = (ico.points[0].x + ico.points[1].x + ico.points[2].x)/3;
    double y2  = (ico.points[0].y + ico.points[1].y + ico.points[2].y)/3;
    double z2  = (ico.points[0].z + ico.points[1].z + ico.points[2].z)/3;
    Dode.points.push_back(Vector3D::point(x2,y2,z2));

    double x3  = (ico.points[0].x + ico.points[3].x + ico.points[4].x)/3;
    double y3  = (ico.points[0].y + ico.points[3].y + ico.points[4].y)/3;
    double z3  = (ico.points[0].z + ico.points[3].z + ico.points[4].z)/3;
    Dode.points.push_back(Vector3D::point(x3,y3,z3));

    double x4  = (ico.points[0].x + ico.points[4].x + ico.points[5].x)/3;
    double y4  = (ico.points[0].y + ico.points[4].y + ico.points[5].y)/3;
    double z4  = (ico.points[0].z + ico.points[4].z + ico.points[5].z)/3;
    Dode.points.push_back(Vector3D::point(x4,y4,z4));

    double x5  = (ico.points[0].x + ico.points[5].x + ico.points[1].x)/3;
    double y5  = (ico.points[0].y + ico.points[5].y + ico.points[1].y)/3;
    double z5  = (ico.points[0].z + ico.points[5].z + ico.points[1].z)/3;
    Dode.points.push_back(Vector3D::point(x5,y5,z5));

    double x6  = (ico.points[1].x + ico.points[6].x + ico.points[2].x)/3;
    double y6  = (ico.points[1].y + ico.points[6].y + ico.points[2].y)/3;
    double z6  = (ico.points[1].z + ico.points[6].z + ico.points[2].z)/3;
    Dode.points.push_back(Vector3D::point(x6,y6,z6));

    double x7  = (ico.points[2].x + ico.points[6].x + ico.points[7].x)/3;
    double y7  = (ico.points[2].y + ico.points[6].y + ico.points[7].y)/3;
    double z7  = (ico.points[2].z + ico.points[6].z + ico.points[7].z)/3;
    Dode.points.push_back(Vector3D::point(x7,y7,z7));

    double x8  = (ico.points[2].x + ico.points[7].x + ico.points[3].x)/3;
    double y8  = (ico.points[2].y + ico.points[7].y + ico.points[3].y)/3;
    double z8  = (ico.points[2].z + ico.points[7].z + ico.points[3].z)/3;
    Dode.points.push_back(Vector3D::point(x8,y8,z8));

    double x9  = (ico.points[3].x + ico.points[7].x + ico.points[8].x)/3;
    double y9  = (ico.points[3].y + ico.points[7].y + ico.points[8].y)/3;
    double z9  = (ico.points[3].z + ico.points[7].z + ico.points[8].z)/3;
    Dode.points.push_back(Vector3D::point(x9,y9,z9));

    double x10  = (ico.points[3].x + ico.points[8].x + ico.points[4].x)/3;
    double y10  = (ico.points[3].y + ico.points[8].y + ico.points[4].y)/3;
    double z10  = (ico.points[3].z + ico.points[8].z + ico.points[4].z)/3;
    Dode.points.push_back(Vector3D::point(x10,y10,z10));


}