#include <fstream>
#include <iostream>
#include <cmath>

#include "easy_image.h"
#include "ini_configuration.h"
#include "LSystems.h"
#include "l_parser/l_parser.h"
#include "vector/vector3d.h"
#include "LSystems.cpp"
#include "ThreeDeeLines.h"
#include "3DFigures.cpp"

#include <string>
#include <list>

#define M_PI 3.14159265358979

img::EasyImage Draw2DLines(Lines2D &lines, const int size, Color bgc) {
    img::Color backgroundcolor(bgc.red*255, bgc.green*255, bgc.blue*255);

    double xmax = 0, xmin = size, ymax = 0, ymin = size;
    for (auto &l: lines){
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
    img::EasyImage image(imagex, imagey, backgroundcolor);

    for (auto& l: lines){
        l.p1.x = (l.p1.x * d) + dcx;
        l.p2.x = (l.p2.x * d) + dcx;
        l.p1.y = (l.p1.y * d) + dcy;
        l.p2.y = (l.p2.y * d) + dcy;
        image.draw_line(
             (unsigned int)std::lround(l.p1.x),
             (unsigned int)std::lround(l.p1.y),
             (unsigned int)std::lround(l.p2.x),
             (unsigned int)std::lround(l.p2.y),
             img::Color(l.color.red*255, l.color.green*255, l.color.blue*255)
        );
    }
    return image;
}



img::EasyImage generate_image(const ini::Configuration &configuration) {
    auto type = configuration["General"]["type"].as_string_or_die();
    auto size = configuration["General"]["size"].as_int_or_die();

    auto backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color bgcolor(backgroundcolor[0], backgroundcolor[1], backgroundcolor[2]);

    if (type == "2DLSystem") {
        auto inputfile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        auto color = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        std::map<int, double> colormap;
        int cm = 0;
        for (auto c: color) {
            colormap[cm] = c;
            cm++;
        }
        Color color2(colormap[0], colormap[1], colormap[2]);

        std::ifstream input(inputfile);
        LParser::LSystem2D lsystem(input);
        input.close();

        auto alfabet = lsystem.get_alphabet();
        auto starting_angle = lsystem.get_starting_angle();
        auto angle = lsystem.get_angle();
        auto initiator = lsystem.get_initiator();
        auto iterations = lsystem.get_nr_iterations();

        double angle_rad = (angle * M_PI) / 180;
        double startingangle_rad = (starting_angle * M_PI) / 180;
        double tempangle;

        int iter = 0;
        std::string full_string = initiator;
        while (iter < iterations) { //Terug veranderen naar iterations!!
            std::string temp;
            for (auto &c: full_string) {
                if (c != '+' && c != '-' && c != '(' && c != ')') {
                    auto replacement = lsystem.get_replacement(c);
                    temp += replacement;
                } else {
                    temp += c;
                }
            }
            full_string = temp;
            iter++;
        }
        //lijnenvector maken
        double x = 0, y = 0;
        Point2D pointend;
        std::list<Point2D> bracket_stack;
        std::list<double> angle_stack;
        Lines2D lijntjes;
        for (auto &c: full_string) {
            if (alfabet.count(c)) {
                if (lsystem.draw(c) == 1) {
                    Point2D p1(x, y);
                    x += std::cos(startingangle_rad);
                    y += std::sin(startingangle_rad);
                    Point2D p2(x, y);
                    Line2D l(p1, p2, color2);
                    lijntjes.push_back(l);
                } else {
                    x += std::cos(startingangle_rad);
                    y += std::sin(startingangle_rad);
                }
            } else if (c == '+') {
                startingangle_rad += angle_rad;
            } else if (c == '-') {
                startingangle_rad -= angle_rad;
            } else if (c == '(') {
                pointend.x = x;
                pointend.y = y;
                angle_stack.push_back(startingangle_rad);
                bracket_stack.push_back(pointend);
            } else if (c == ')') {
                Point2D last_point = bracket_stack.back();
                startingangle_rad = angle_stack.back();
                x = last_point.x;
                y = last_point.y;
                bracket_stack.pop_back();
                angle_stack.pop_back();
            }
        }
        auto image = Draw2DLines(lijntjes, size, bgcolor);
        return image;
    }

    else if (type == "Wireframe") {
        Figures3D figures;
        auto nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

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
            else if (type2 == "Cube"){
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
            else if (type2 == "Tetrahedron"){
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
            else if (type2 == "Octahedron"){
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
            else if (type2 == "Icosahedron"){
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
            else if (type2 == "Dodecahedron"){
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
            else if (type2 == "Sphere"){
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
            else if (type2 == "Cone"){
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
            else if (type2 == "Cylinder"){
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
            else if (type2 == "Torus"){
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

                auto sphere = generateTorus(figureColor, r, R, n ,m);
                applyTransformation(sphere, allTrans);
                figures.push_back(sphere);
            }
            else if (type2 == "3DLSystem"){
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

                angle = angle*M_PI/180;

                Vector3D curp = Vector3D::point(0,0,0);
                Vector3D H = Vector3D::point(1,0,0);
                Vector3D L = Vector3D::point(0,1,0);
                Vector3D U = Vector3D::point(0,0,1);

                int counter = 0;
                std::vector<Vector3D> pointstack;
                std::vector<Vector3D> Hstack;
                std::vector<Vector3D> Lstack;
                std::vector<Vector3D> Ustack;

                fig.points.push_back(curp);

                for (int i = 0; i < iterations; i++){
                    std::string temp;
                    for (auto s: initiator){
                        if (s == '+') temp+="+";
                        else if (s == '-') temp+="-";
                        else if (s == '^') temp +="^";
                        else if (s == '&') temp+="&";
                        else if (s == '\\') temp +="\\"
                    }
                }
            }

        }
        auto lijntjes = doProjection(figures);
        auto img = Draw2DLines(lijntjes, size, bgcolor);
        return img;
    }
    return img::EasyImage();
}


int main(int argc, char const *argv[]) {
    int retVal = 0;
    try
    {
        std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for(std::string fileName : args)
        {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(fileName);
                if (fin.peek() == std::istream::traits_type::eof()) {
                    std::cout << "Ini file appears empty. Does '" <<
                              fileName << "' exist?" << std::endl;
                    continue;
                }
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}
