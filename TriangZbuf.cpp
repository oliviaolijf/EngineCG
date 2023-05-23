#include "ThreeDeeLines.h"
#include "Zbuffer.h"
#include "LSystems.h"
#include "3DFigures.cpp"
#include "ini_configuration.h"
#include "l_parser/l_parser.h"
#include <fstream>
#include <limits>

std::vector<Face> triangulate(const Face& face){
    auto n = face.point_indexes.size();
    std::vector<Face> triangles;

    for (int i = 2; i < n-2; i++){
        Face f = {{face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]}};
        triangles.push_back(f);
    }
    return triangles;
}

void draw_zbuf_triag(ZBuffer& buffer, img::EasyImage& img,
                     Vector3D const& A,
                     Vector3D const& B,
                     Vector3D const& C,
                     double d, double dx, double dy,
                     Color color){
    img::Color pixelcolor(color.red*255, color.green*255, color.blue*255);
    Point2D aAccent, bAccent, cAccent;
    auto xa = (d*A.x/-A.z) + dx;
    auto ya = (d*A.y/-A.z) + dy;

    auto xb = (d*B.x/-B.z) + dx;
    auto yb = (d*B.y/-B.z) + dy;

    auto xc = (d*C.x/-C.z) + dx;
    auto yc = (d*C.y/-C.z) + dy;

    int ymin = std::min(ya, yb);
    if (yc < ymin) ymin = yc;
    int ymax = std::max(ya, yb);
    if (yc > ymax) ymax = yc;
    int xmin = std::min(xa, xb);
    if (xc < xmin) xmin = xc;
    int xmax = std::max(xa, xb);
    if (xc > xmax) xmax = xc;

    for (auto i = 0; i <= (ymax-ymin); i++){
        double xlAB = std::numeric_limits<double>::infinity();
        double xlAC = std::numeric_limits<double>::infinity();
        double xlBC = std::numeric_limits<double>::infinity();
        double xrAB = -std::numeric_limits<double>::infinity();
        double xrAC = -std::numeric_limits<double>::infinity();
        double xrBC = -std::numeric_limits<double>::infinity();
        int yi = ymin+i;

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
        double templ = round(std::min(xlAB, xlAC));
        int xl = round(std::min(templ, xlBC) + 0.5);

        double tempr = round(std::min(xrAB, xrAC));
        int xr = round(std::min(tempr, xrBC) - 0.5);

        for (int ix = xl; ix <= xr; ix++ ){
            img(ix, yi) = pixelcolor;
        }
    }
}

img::EasyImage triangZbuf(const ini::Configuration& configuration){
    Figures3D figures;
    Figures3D triFigures;
    auto nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    auto size = configuration["General"]["size"].as_int_or_die();

    auto backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color bgcolor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);

    std::vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eyepoint = Vector3D::point(eye[0], eye[1], eye[2]);
    auto V = eyePointTrans(eyepoint);

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
            fig.color = drawingcolor;
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto cube = generateCube(figureColor);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto cube = generateTetrahedron(figureColor);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto cube = generateOctahedron(figureColor);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto cube = generateIcosahedron(figureColor);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto cube = generateDodecahedron(figureColor);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto sphere = generateSphere(figureColor, iterations);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto sphere = generateCone(figureColor, iterations, height);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto sphere = generateCylinder(figureColor, iterations, height);
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

            auto figure_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            Color figureColor = Color(figure_color[0], figure_color[1], figure_color[2]);

            auto sphere = generateTorus(figureColor, r, R, n, m);
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
                    keepp =curp;
                    keepH = H;
                    keepL = L;
                    keepU = U;
                    /**pointstack.push_back(curp);
                    Hstack.push_back(H);
                    Lstack.push_back(L);
                    Ustack.push_back(U);**/
                    continue;
                }
                if (c == ')'){
                    H = keepH;
                    L = keepL;
                    U = keepU;
                    curp = keepp;
                    /**curp = pointstack.back();
                    pointstack.pop_back();

                    H =  Hstack.back();
                    Hstack.pop_back();

                    L = Lstack.back();
                    Lstack.pop_back();

                    U = Ustack.back();
                    Ustack.pop_back();**/
                    continue;
                }
            }
            fig.color = figureColor;
            applyTransformation(fig, allTrans);
            figures.push_back(fig);
        }
    }

    for (auto fig: figures){
        Figure temp;
        for (auto &face: fig.faces){
            auto huh = triangulate(face);
            for (auto h: huh)  temp.faces.push_back(h);
        }
        triFigures.push_back(temp);
    }
    auto lijntjes = doProjection(triFigures);

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

    auto DCx = (d*(xmin+xmax))/2;
    auto DCy = (d*(ymin+ymax))/2;
    auto dcx = (imagex/2)-DCx;
    auto dcy = (imagey/2)-DCy;


}