#include "ThreeDeeLines.h"
#include "Zbuffer.h"
#include "ini_configuration.h"
#include "3DFigures.cpp"
#include "l_parser/l_parser.h"
#include <fstream>
#include <limits>

#define M_PI 3.14159265358979

std::vector<Face> triangulate(const Face& face){
    auto n = face.point_indexes.size();
    std::vector<Face> triangles;

    if (face.point_indexes.size() > 3) {
        for (int i = 1; i < n - 1; i++) {
            Face f = {{face.point_indexes[0], face.point_indexes[i], face.point_indexes[i + 1]}};
            triangles.push_back(f);
        }
    }
    else triangles.push_back(face);
    return triangles;
}

void draw_zbuf_triag(ZBuffer& buffer, img::EasyImage& img,
                     Vector3D const& A,
                     Vector3D const& B,
                     Vector3D const& C,
                     double d, double dx, double dy,
                     Color ambientReflection,
                     Color diffuseReflection,
                     Color specularReflection,
                     Light light, double reflectionCoeff)
                     {
    Color am (light.ambientLight.red * ambientReflection.red,
                light.ambientLight.green * ambientReflection.green,
                light.ambientLight.blue * ambientReflection.blue
            );
    img::Color pixelcolor (am.red*255, am.green*255, am.blue*255);

    auto xa = ((d*A.x)/-A.z) + dx;
    auto ya = ((d*A.y)/-A.z) + dy;
    auto xb = ((d*B.x)/-B.z) + dx;
    auto yb = ((d*B.y)/-B.z) + dy;
    auto xc = ((d*C.x)/-C.z) + dx;
    auto yc = ((d*C.y)/-C.z) + dy;

    double tempymin = std::min(ya, yb);
    int ymin = std::lround(std::min(tempymin, yc)+0.5);
    double tempymax = std::max(ya, yb);
    int ymax = std::lround(std::max(tempymax, yc)-0.5);

    for (auto yi = ymin; yi <= ymax; yi++){
        double xlAB = std::numeric_limits<double>::infinity();
        double xlAC = std::numeric_limits<double>::infinity();
        double xlBC = std::numeric_limits<double>::infinity();
        double xrAB = -std::numeric_limits<double>::infinity();
        double xrAC = -std::numeric_limits<double>::infinity();
        double xrBC = -std::numeric_limits<double>::infinity();

        double xp = xa;
        double xq = xb;
        double yp = ya;
        double yq = yb;

        if ((yi-yp)*(yi-yq)<= 0 && yp != yq) {
            double xi = xq + (xp - xq) * (yi - yq) / (yp - yq);
            xlAB = xi;
            xrAB = xi;
        }

        xq = xc;
        yq = yc;
        if ((yi-yp)*(yi-yq)<= 0 && yp != yq) {
            double xi = xq + (xp - xq) * (yi - yq) / (yp - yq);
            xlAC = xi;
            xrAC = xi;
        }

        xp = xb;
        yp = yb;
        if ((yi-yp)*(yi-yq)<= 0 && yp != yq) {
            double xi = xq + (xp - xq) * (yi - yq) / (yp - yq);
            xlBC = xi;
            xrBC = xi;
        }
        double templ = std::min(xlAB, xlAC);
        int xl = std::lround(std::min(templ, xlBC) + 0.5);
        double tempr = std::max(xrAB, xrAC);
        int xr = std::lround(std::max(tempr, xrBC) - 0.5);

        auto xg = (xa+xb+xc)/3;
        auto yg = (ya+yb+yc)/3;
        double smallzg = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);

        Vector3D u = B-A;
        Vector3D v = C-A;
        double w1 = (u.y*v.z) - (u.z*v.y);
        double w2 = (u.z*v.x) - (u.x*v.z);
        double w3 = (u.x*v.y) - (u.y*v.x);
        double k = w1*A.x + w2*A.y + w3*A.z;

        double dzdx = w1/(-d*k);
        double dzdy = w2/(-d*k);

        for (int ix = xl; ix <= xr; ix++ ){
            auto curz = buffer[ix][yi];
            auto newz = 1.0001 * smallzg + (ix-xg)*dzdx + (yi-yg)*dzdy;
            if (newz < curz) {
                img(ix, yi) = pixelcolor;
                buffer[ix][yi] = newz;
            }
        }
    }
}

img::EasyImage triangZbuf(const ini::Configuration& configuration){
    Figures3D figures;
    Figures3D triFigures;
    auto type = configuration["General"]["type"].as_string_or_die();
    auto nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    auto size = configuration["General"]["size"].as_int_or_die();

    auto backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    img::Color bgcolor(backgroundcolor[0]*255, backgroundcolor[1]*255, backgroundcolor[2]*255);

    std::vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eyepoint = Vector3D::point(eye[0], eye[1], eye[2]);
    auto V = eyePointTrans(eyepoint);
    Color white(1,1,1);
    Light finalLight;

    if (type == "LightedZBuffering"){
        auto nrl = configuration["General"]["nrLights"].as_int_or_die();
        Lights3D lights;

        ini::DoubleTuple defa = {0, 0, 0};
        for (int i = 0; i < nrl; i++){
            Light light;
            std::string lightname = "Light"+ std::to_string(i);

            auto al = configuration[lightname]["ambientLight"].as_double_tuple_or_default(defa);
            Color ambient(al[0], al[1], al[2]);
            light.ambientLight = ambient;

            auto dl = configuration[lightname]["diffuseLight"].as_double_tuple_if_exists(defa);
            Color diffuse(al[0], al[1], al[2]);
            light.diffuseLight = diffuse;

            auto sl = configuration[lightname]["specularLight"].as_double_tuple_if_exists(defa);
            Color specular(al[0], al[1], al[2]);
            light.specularLight = specular;

            lights.push_back(light);
        }

        for (auto l: lights){
            finalLight.ambientLight.red += l.ambientLight.red;
            finalLight.ambientLight.green += l.ambientLight.green;
            finalLight.ambientLight.blue += l.ambientLight.blue;

            finalLight.diffuseLight.red += l.diffuseLight.red;
            finalLight.diffuseLight.green += l.diffuseLight.green;
            finalLight.diffuseLight.blue += l.diffuseLight.blue;

            finalLight.specularLight.red += l.specularLight.red;
            finalLight.specularLight.green += l.specularLight.green;
            finalLight.specularLight.blue += l.specularLight.blue;
        }

    }

    for (int i = 0; i < nrFigures; i++) {
        std::string figure_name = "Figure" + std::to_string(i);
        auto type2 = configuration[figure_name]["type"].as_string_or_die();

        if (type2 == "LineDrawing") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto nrpoints = configuration[figure_name]["nrPoints"].as_int_or_die();
            auto nrlines = configuration[figure_name]["nrLines"].as_int_or_die();

            Figure fig;
            for (int j = 0; j < nrpoints; j++) {
                std::string point_name = "point" + std::to_string(j);
                std::vector<double> point = configuration[figure_name][point_name].as_double_tuple_or_die();
                fig.points.push_back(Vector3D::point(point[0], point[1], point[2]));
            }
            for (int j = 0; j < nrlines; j++) {
                std::string line_name = "line" + std::to_string(j);
                std::vector<int> line = configuration[figure_name][line_name].as_int_tuple_or_die();
                Face face;
                for (auto a: line) { face.point_indexes.push_back(a); }
                fig.faces.push_back(face);
            }
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                fig.ambientReflection = ambiR;

                fig.specularReflection = Color(0,0,0);
                fig.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = drawingcolor;
                fig.ambientReflection = Color(1,1,1);
                fig.specularReflection = Color(0,0,0);
                fig.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(fig, allTrans);
            figures.push_back(fig);
        }
        else if (type2 == "Cube") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;


            auto cube = generateCube(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(cube, allTrans);
            figures.push_back(cube);
        }
        else if (type2 == "Tetrahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateTetrahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(cube, allTrans);
            figures.push_back(cube);
        }
        else if (type2 == "Octahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateOctahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }

            applyTransformation(cube, allTrans);
            figures.push_back(cube);
        }
        else if (type2 == "Icosahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateIcosahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(cube, allTrans);
            figures.push_back(cube);
        }
        else if (type2 == "Dodecahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateDodecahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(cube, allTrans);
            figures.push_back(cube);
        }
        else if (type2 == "Sphere") {
            auto iterations = configuration[figure_name]["n"].as_int_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateSphere(white, iterations);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color (1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(sphere, allTrans);
            figures.push_back(sphere);
        }
        else if (type2 == "Cone") {
            auto iterations = configuration[figure_name]["n"].as_int_or_die();
            auto height = configuration[figure_name]["height"].as_double_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateCone(white, iterations, height);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(sphere, allTrans);
            figures.push_back(sphere);
        }
        else if (type2 == "Cylinder") {
            auto iterations = configuration[figure_name]["n"].as_int_or_die();
            auto height = configuration[figure_name]["height"].as_double_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateCylinder(white, iterations, height);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(sphere, allTrans);
            figures.push_back(sphere);
        }
        else if (type2 == "Torus") {
            auto r = configuration[figure_name]["r"].as_double_or_die();
            auto R = configuration[figure_name]["R"].as_double_or_die();
            auto n = configuration[figure_name]["n"].as_double_or_die();
            auto m = configuration[figure_name]["m"].as_double_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateTorus(white, r, R, n, m);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(sphere, allTrans);
            figures.push_back(sphere);
        }
        else if (type2 == "3DLSystem") {
            Figure fig;
            auto inputfile = configuration[figure_name]["inputfile"].as_string_or_die();
            std::ifstream input(inputfile);
            LParser::LSystem3D lsystem(input);
            input.close();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto alfabet = lsystem.get_alphabet();
            std::string initiator = lsystem.get_initiator();
            auto angle = lsystem.get_angle();
            auto iterations = lsystem.get_nr_iterations();

            angle = (angle * M_PI) / 180;

            int pointcounter = 0;

            Vector3D curp = Vector3D::point(0, 0, 0);
            Vector3D H = Vector3D::point(1, 0, 0);
            Vector3D L = Vector3D::point(0, 1, 0);
            Vector3D U = Vector3D::point(0, 0, 1);

            Vector3D keepH;
            Vector3D keepL;
            Vector3D keepU;
            Vector3D keepp;

            fig.points.push_back(curp);

            int counter = 0;
            std::list<Vector3D> pointstack;
            std::list<Vector3D> Hstack;
            std::list<Vector3D> Lstack;
            std::list<Vector3D> Ustack;

            for (int it = 0; it < iterations; it++) {
                std::string temp;
                for (auto &s: initiator) {
                    if (s == '+') temp += "+";
                    else if (s == '-') temp += "-";
                    else if (s == '^') temp += "^";
                    else if (s == '&') temp += "&";
                    else if (s == '\\') temp += "\\";
                    else if (s == '/') temp += "/";
                    else if (s == '|') temp += "|";
                    else if (alfabet.count(s)) temp += lsystem.get_replacement(s);
                }
                initiator = temp;
            }
            for (auto &c: initiator) {
                if (alfabet.count(c)) {
                    if (lsystem.draw(c)) {
                        fig.points.push_back(curp);
                        curp += H;
                        pointcounter++;
                        Face face({pointcounter, pointcounter-1});
                        fig.faces.push_back(face);
                        fig.points.push_back(curp);
                    } else {
                        curp += H;
                    }
                    continue;
                }
                if (c == '+') {
                    Vector3D temp = H;
                    H = (H*cos(angle)) + (L*sin(angle));
                    L = -(temp*sin(angle)) + (L*cos(angle));
                    continue;
                }
                if (c == '-') {
                    Vector3D temp = H;
                    H = (H*cos(-angle)) + (L* sin(-angle));
                    L = -(temp*sin(-angle)) + (L*cos(-angle));
                    continue;
                }
                if (c == '^'){
                    Vector3D temp = H;
                    H = (H*cos(angle)) + (U* sin(angle));
                    U = -(temp*sin(angle)) + (U*cos(angle));
                    continue;
                }
                if (c == '&') {
                    Vector3D temp = H;
                    H = (H* cos(-angle)) + (U* sin(-angle));
                    U = -(temp* sin(-angle)) + (U*cos(-angle));
                    continue;
                }
                if (c == '\\'){
                    Vector3D temp = L;
                    L = (L* cos(angle)) - (U*sin(angle));
                    U = (temp*sin(angle)) + (U* cos(angle));
                    continue;
                }
                if (c == '/'){
                    Vector3D temp = L;
                    L = (L*cos(-angle)) - (U* sin(-angle));
                    U = (temp*sin(-angle)) + (U*cos(-angle));
                    continue;
                }
                if (c == '|'){
                    H = -H;
                    L = -L;
                    continue;
                }
                if (c == '('){
                    pointstack.push_back(curp);
                    Hstack.push_back(H);
                    Lstack.push_back(L);
                    Ustack.push_back(U);
                    continue;
                }
                if (c == ')'){
                    curp = pointstack.back();
                    pointstack.pop_back();

                    H =  Hstack.back();
                    Hstack.pop_back();

                    L = Lstack.back();
                    Lstack.pop_back();

                    U = Ustack.back();
                    Ustack.pop_back();
                    continue;
                }
            }
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                fig.ambientReflection = ambiR;

                fig.specularReflection = Color(0,0,0);
                fig.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                fig.ambientReflection = drawingcolor;
                fig.specularReflection = Color(0,0,0);
                fig.diffuseReflection = Color(0,0,0);
            }
            applyTransformation(fig, allTrans);
            figures.push_back(fig);
        }
        else if (type2 == "ThickLineDrawing") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto nrpoints = configuration[figure_name]["nrPoints"].as_int_or_die();
            auto nrlines = configuration[figure_name]["nrLines"].as_int_or_die();

            std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);

            Figure fig;
            for (int j = 0; j < nrpoints; j++) {
                std::string point_name = "point" + std::to_string(j);
                std::vector<double> point = configuration[figure_name][point_name].as_double_tuple_or_die();
                fig.points.push_back(Vector3D::point(point[0], point[1], point[2]));
            }
            for (int j = 0; j < nrlines; j++) {
                std::string line_name = "line" + std::to_string(j);
                std::vector<int> line = configuration[figure_name][line_name].as_int_tuple_or_die();
                Face face;
                for (auto a: line) { face.point_indexes.push_back(a); }
                fig.faces.push_back(face);
            }
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                fig.ambientReflection = ambiR;

                fig.specularReflection = Color(0,0,0);
                fig.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                fig.ambientReflection = drawingcolor;
                fig.specularReflection = Color(0,0,0);
                fig.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(fig, thickfigs, radius, n, m);
            for (auto &f: thickfigs){
                applyTransformation(f, allTrans);
                figures.push_back(f);
            }
        }
        else if (type2 == "ThickCube") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateCube(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(cube, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickTetrahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto tetra = generateTetrahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                tetra.ambientReflection = ambiR;

                tetra.specularReflection = Color(0,0,0);
                tetra.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                tetra.ambientReflection = drawingcolor;
                tetra.specularReflection = Color(0,0,0);
                tetra.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(tetra, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickOctahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateOctahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(cube, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickIcosahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateIcosahedron(white);
            if (type == "LightedZBuffering") {
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;
                cube.specularReflection = Color(0, 0, 0);
                cube.diffuseReflection = Color(0, 0, 0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(cube, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickDodecahedron") {
            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();

            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto cube = generateDodecahedron(white);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                cube.ambientReflection = ambiR;

                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                cube.ambientReflection = drawingcolor;
                cube.specularReflection = Color(0,0,0);
                cube.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(cube, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickSphere") {
            auto iterations = configuration[figure_name]["n"].as_int_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();
            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateSphere(white, iterations);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(sphere, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickCone") {
            auto iterations = configuration[figure_name]["n"].as_int_or_die();
            auto height = configuration[figure_name]["height"].as_double_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();
            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateCone(white, iterations, height);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(sphere, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickCylinder") {
            auto iterations = configuration[figure_name]["n"].as_int_or_die();
            auto height = configuration[figure_name]["height"].as_double_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();
            auto m = configuration[figure_name]["m"].as_int_or_die();
            auto n = configuration[figure_name]["n"].as_int_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateCylinder(white, iterations, height);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(sphere, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
        else if (type2 == "ThickTorus") {
            auto r = configuration[figure_name]["r"].as_double_or_die();
            auto R = configuration[figure_name]["R"].as_double_or_die();
            auto n = configuration[figure_name]["n"].as_double_or_die();
            auto m = configuration[figure_name]["m"].as_double_or_die();

            auto rotatX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotatY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotatZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scal = configuration[figure_name]["scale"].as_double_or_die();
            auto radius = configuration[figure_name]["radius"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto centerr = Vector3D::point(center[0], center[1], center[2]);
            auto trans = translate(centerr);

            auto rotX = rotateX((rotatX * M_PI) / 180);
            auto rotY = rotateY((rotatY * M_PI) / 180);
            auto rotZ = rotateZ((rotatZ * M_PI) / 180);
            auto scale = scaleFigure(scal);
            auto allTrans = rotX * rotZ * rotY * scale * trans * V;

            auto sphere = generateTorus(white, r, R, n, m);
            if (type == "LightedZBuffering"){
                auto ambi_R = configuration[figure_name]["ambientReflection"].as_double_tuple_or_die();
                Color ambiR(ambi_R[0], ambi_R[1], ambi_R[2]);
                sphere.ambientReflection = ambiR;

                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            else {
                std::vector<double> drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
                Color drawingcolor = Color(drawing_color[0], drawing_color[1], drawing_color[2]);
                finalLight.ambientLight = Color(1,1,1);
                sphere.ambientReflection = drawingcolor;
                sphere.specularReflection = Color(0,0,0);
                sphere.diffuseReflection = Color(0,0,0);
            }
            Figures3D thickfigs;
            generateThickFigure(sphere, thickfigs, radius, n, m);
            for (auto &fig: thickfigs){
                applyTransformation(fig, allTrans);
                figures.push_back(fig);
            }
        }
    }
    for (auto& fig: figures){
        Figure temp;
        for (auto &face: fig.faces){
            auto huh = triangulate(face);
            for (auto h: huh)  temp.faces.push_back(h);
        }
        fig.faces = temp.faces;
    }
    auto copyfigures = figures;
    auto lijntjes = doProjection(copyfigures);

    double xmax = 0, xmin = size, ymax = 0, ymin = size;
    for (auto &l: lijntjes){
        if (std::max(l.p1.x, l.p2.x) > xmax) xmax = std::max(l.p2.x, l.p1.x);
        if (std::min(l.p1.x, l.p2.x) < xmin) xmin = std::min(l.p2.x, l.p1.x);
        if (std::max(l.p1.y, l.p2.y) > ymax) ymax = std::max(l.p2.y, l.p1.y);
        if (std::min(l.p1.y, l.p2.y) < ymin) ymin = std::min(l.p2.y, l.p1.y);
    }

    auto xrange = xmax-xmin;
    auto yrange = ymax-ymin;

    auto imagex = size*(xrange/ std::max(xrange, yrange));
    auto imagey = size*(yrange/ std::max(xrange, yrange));

    auto d = 0.95*(imagex/xrange);

    auto DCx = d*((xmin+xmax)/2);
    auto DCy = d*((ymin+ymax)/2);
    auto dcx = (imagex/2)-DCx;
    auto dcy = (imagey/2)-DCy;

    img::EasyImage image(std::lround(imagex), std::lround(imagey), bgcolor);

    ZBuffer buffer(image.get_width(), image.get_height());
    double idk = 1;
    for (auto& fig: figures){
        for (auto& fa: fig.faces){
                Vector3D a = fig.points[fa.point_indexes[0]];
                Vector3D b = fig.points[fa.point_indexes[1]];
                Vector3D c = fig.points[fa.point_indexes[2]];

                draw_zbuf_triag(buffer, image, a, b, c, d, dcx, dcy, fig.ambientReflection, fig.diffuseReflection,
                                fig.specularReflection, finalLight, idk);
        }
    }
    return image;
}