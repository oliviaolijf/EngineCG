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

    double x2  = (ico.points[0].x + ico.points[2].x + ico.points[3].x)/3;
    double y2  = (ico.points[0].y + ico.points[2].y + ico.points[3].y)/3;
    double z2  = (ico.points[0].z + ico.points[2].z + ico.points[3].z)/3;
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

    double x10 = (ico.points[3].x + ico.points[8].x + ico.points[4].x)/3;
    double y10 = (ico.points[3].y + ico.points[8].y + ico.points[4].y)/3;
    double z10 = (ico.points[3].z + ico.points[8].z + ico.points[4].z)/3;
    Dode.points.push_back(Vector3D::point(x10,y10,z10));

    double x11 = (ico.points[4].x + ico.points[8].x + ico.points[9].x)/3;
    double y11 = (ico.points[4].y + ico.points[8].y + ico.points[9].y)/3;
    double z11 = (ico.points[4].z + ico.points[8].z + ico.points[9].z)/3;
    Dode.points.push_back(Vector3D::point(x11,y11,z11));

    double x12 = (ico.points[4].x + ico.points[9].x + ico.points[5].x)/3;
    double y12 = (ico.points[4].y + ico.points[9].y + ico.points[5].y)/3;
    double z12 = (ico.points[4].z + ico.points[9].z + ico.points[5].z)/3;
    Dode.points.push_back(Vector3D::point(x12,y12,z12));

    double x13 = (ico.points[5].x + ico.points[9].x + ico.points[10].x)/3;
    double y13 = (ico.points[5].y + ico.points[9].y + ico.points[10].y)/3;
    double z13 = (ico.points[5].z + ico.points[9].z + ico.points[10].z)/3;
    Dode.points.push_back(Vector3D::point(x13,y13,z13));

    double x14 = (ico.points[5].x + ico.points[10].x + ico.points[1].x)/3;
    double y14 = (ico.points[5].y + ico.points[10].y + ico.points[1].y)/3;
    double z14 = (ico.points[5].z + ico.points[10].z + ico.points[1].z)/3;
    Dode.points.push_back(Vector3D::point(x14,y14,z14));

    double x15 = (ico.points[1].x + ico.points[10].x + ico.points[6].x)/3;
    double y15 = (ico.points[1].y + ico.points[10].y + ico.points[6].y)/3;
    double z15 = (ico.points[1].z + ico.points[10].z + ico.points[6].z)/3;
    Dode.points.push_back(Vector3D::point(x15,y15,z15));

    double x16 = (ico.points[11].x + ico.points[7].x + ico.points[6].x)/3;
    double y16 = (ico.points[11].y + ico.points[7].y + ico.points[6].y)/3;
    double z16 = (ico.points[11].z + ico.points[7].z + ico.points[6].z)/3;
    Dode.points.push_back(Vector3D::point(x16,y16,z16));

    double x17 = (ico.points[11].x + ico.points[8].x + ico.points[7].x)/3;
    double y17 = (ico.points[11].y + ico.points[8].y + ico.points[7].y)/3;
    double z17 = (ico.points[11].z + ico.points[8].z + ico.points[7].z)/3;
    Dode.points.push_back(Vector3D::point(x17,y17,z17));

    double x18 = (ico.points[11].x + ico.points[9].x + ico.points[8].x)/3;
    double y18 = (ico.points[11].y + ico.points[9].y + ico.points[8].y)/3;
    double z18 = (ico.points[11].z + ico.points[9].z + ico.points[8].z)/3;
    Dode.points.push_back(Vector3D::point(x18,y18,z18));

    double x19 = (ico.points[11].x + ico.points[10].x + ico.points[9].x)/3;
    double y19 = (ico.points[11].y + ico.points[10].y + ico.points[9].y)/3;
    double z19 = (ico.points[11].z + ico.points[10].z + ico.points[9].z)/3;
    Dode.points.push_back(Vector3D::point(x19,y19,z19));

    double x20 = (ico.points[11].x + ico.points[6].x + ico.points[10].x)/3;
    double y20 = (ico.points[11].y + ico.points[6].y + ico.points[10].y)/3;
    double z20 = (ico.points[11].z + ico.points[6].z + ico.points[10].z)/3;
    Dode.points.push_back(Vector3D::point(x20,y20,z20));

    Face f1({0, 1, 2, 3, 4});
    Face f2({0, 5, 6, 7, 1});
    Face f3({1, 7, 8, 9, 2});
    Face f4({2, 9, 10, 11, 3});
    Face f5({3, 11, 12, 13, 4});
    Face f6({4, 13, 14, 5, 0});
    Face f7({19, 18, 17, 16, 15});
    Face f8({19, 14, 13, 12, 18});
    Face f9({18, 12, 11, 10, 17});
    Face f10({17, 10, 9, 8, 16});
    Face f11({16, 8, 7, 6, 15});
    Face f12({15, 6, 5, 14, 19});

    Dode.faces = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12};

    Dode.color = c;
    return Dode;
}

Figure generateSphere(Color c, const int n){
    auto ico = generateIcosahedron(c);
    int points = 0;
    Figure sphere;

    for (int i = 0; i < n; i++) {
        for (auto &face: ico.faces) {
            Vector3D A = ico.points[face.point_indexes[0]];
            Vector3D B = ico.points[face.point_indexes[1]];
            Vector3D C = ico.points[face.point_indexes[2]];

            Vector3D D = (A + B) / 2;
            Vector3D E = (A + C) / 2;
            Vector3D F = (B + C) / 2;

            sphere.points.push_back(A);
            sphere.points.push_back(B);
            sphere.points.push_back(C);
            sphere.points.push_back(D);
            sphere.points.push_back(E);
            sphere.points.push_back(F);

            Face f1({points, points + 3, points + 4});
            Face f2({points + 1, points + 5, points + 3});
            Face f3({points + 2, points + 4, points + 5});
            Face f4({points + 3, points + 5, points + 4});

            sphere.faces.push_back(f1);
            sphere.faces.push_back(f2);
            sphere.faces.push_back(f3);
            sphere.faces.push_back(f4);

            points += 6;
        }
        ico = sphere;
    }
    for (auto &p: ico.points){
        auto r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        Vector3D newp = Vector3D::point(p.x/r, p.y/r, p.z/r);
        p = newp;
        p.normalise();
    }
    ico.color = c;
    return ico;
}

Figure generateCone(Color c, int n, double height){
    Figure cone;
    for (int i = 0; i < n; i++){
        Vector3D p = Vector3D::point(std::cos(2*i*M_PI/n), std::sin(2*i*M_PI/n), 0);
        cone.points.push_back(p);
    }
    Vector3D p = Vector3D::point(0,0,height);
    cone.points.push_back(p);
    Face last;
    for (int i = 0; i<n ; i++){
        Face f({i, (i+1)%n, n});
        cone.faces.push_back(f);

        last.point_indexes.push_back(n-i-1);
    }
    cone.faces.push_back(last);
    cone.color = c;
    return cone;
}

Figure generateCylinder(Color c, int n, double height){
    Figure cyl;

    for (int i = 0; i < n; i++){
        Vector3D p = Vector3D::point(std::cos(2*i*M_PI/n), std::sin(2*i*M_PI/n), 0);
        cyl.points.push_back(p);
    }
    for (int i = 0; i < n; i++){
        Vector3D p = Vector3D::point(std::cos(2*i*M_PI/n), std::sin(2*i*M_PI/n), height);
        cyl.points.push_back(p);
    }

    for (int i = 0; i < n-1; i++){
        Face f({i, (i+1)%n, (i+n+1)%(2*n), i+n});
        cyl.faces.push_back(f);
    }
    Face laatste({n-1, 0, n, (2*n)-1});
    cyl.faces.push_back(laatste);

    Face onder;
    for (auto i = 0; i<n; i++){onder.point_indexes.push_back(i);}
    cyl.faces.push_back(onder);

    Face boven;
    for (auto i = n; i<2*n; i++){boven.point_indexes.push_back(i);}
    cyl.faces.push_back(boven);

    cyl.color = c;
    return cyl;
}

Figure generateTorus(Color c, double r, double R, int n, int m){
    Figure torus;

    std::vector<std::vector<int>> pointmap(n, std::vector<int>(m,0));
    int points = 0;

    for (int i = 0; i<n; i++){
        for (int j = 0; j<m; j++){
            auto u = (2*i*M_PI)/n;
            auto v = (2*j*M_PI)/m;

            auto x = (R + r*std::cos(v))*std::cos(u);
            auto y = (R + r*std::cos(v))*std::sin(u);
            auto z = r*std::sin(v);

            Vector3D p = Vector3D::point(x,y,z);
            torus.points.push_back(p);

            pointmap[i][j] = points;
            points++;
        }
    }

    for (int i = 0; i<n; i++){
        for (auto j = 0; j<m; j++){
            Face f({pointmap[i][j], pointmap[(i+1)%n][j], pointmap[(i + 1)%n][(j+1)%m], pointmap[i][(j+1)%m]});
            torus.faces.push_back(f);
        }
    }

    torus.color = c;
    return torus;
}

void generateThickFigure(const Figure& fig, Figures3D &result, const double r, const int n, const int m){
    for (auto &face: fig.faces) {
        for (auto &i: face.point_indexes) {
            auto sphere = generateSphere(fig.color, m);
            Vector3D p = Vector3D::vector(fig.points[i].x, fig.points[i].y, fig.points[i].z);
            auto S = scaleFigure(r);
            auto T = translate(p);
            auto ST = S * T;
            applyTransformation(sphere, ST);
            result.push_back(sphere);
        }
        for (int i = 0; i < face.point_indexes.size(); i++){
            Vector3D p1;
            Vector3D p2;
            if (i == face.point_indexes.size() - 1){
                p1 = fig.points[face.point_indexes[i]];
                p2 = fig.points[face.point_indexes[0]];
            }
            else {
                p1 = fig.points[face.point_indexes[i]];
                p2 = fig.points[face.point_indexes[i+1]];
            }
            Vector3D vec(p2-p1);
            double l = vec.length();
            double h = l/r;
            auto cylinder = generateCylinder(fig.color, n, h);

            auto S = scaleFigure(r);
            auto Pr = Vector3D::point(0,0,0) + vec;
            double phi, theta, rr;
            toPolar(Pr, theta, phi, rr);

            auto P = rotateY(phi);
            auto T = rotateZ(theta);
            auto TRA = translate(p1);
            auto alltrans = S*P*T*TRA;

            applyTransformation(cylinder, alltrans);
            result.push_back(cylinder);
        }
    }
}

void generateFractal(Figure& fig, Figures3D& fractal, const int nrit, const double scale){
    Figures3D toconvert;
    toconvert.push_back(fig);
    Figures3D newfigs;
    if (nrit != 0) {
        for (unsigned int it = 0; it < nrit; it++) {
            newfigs.clear();
            for (auto &f: toconvert) {
                for (unsigned int i = 0; i < f.points.size(); i++) {
                    auto curpoint = f.points[i];

                    auto newfig = f;
                    auto S = scaleFigure(1 / scale);
                    applyTransformation(newfig, S);

                    auto T = translate(curpoint - newfig.points[i]);
                    applyTransformation(newfig, T);

                    newfigs.push_back(newfig);
                }
            }
            toconvert = newfigs;
        }
        fractal = newfigs;
    }
    else{
        fractal.push_back(fig);
    }
}