#include "easy_image.h"
#include "ini_configuration.h"
#include "LSystems.h"
#include "l_parser/l_parser.h"
#include "vector/vector3d.h"

#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <string>
#include <list>
#include <cmath>

using Lines2D = std::list<Line2D>;

img::EasyImage Draw2DLines(Lines2D &lines, const int size, Color bgc, Color drawing_color) {
    auto backgroundcolor = bgc.createColor();
    auto drawingcolor = drawing_color.createColor();

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
             drawingcolor
        );
    }
    return image;
}

Point2D doProjectionPoint(const Vector3D &point, const double &d){
    auto xAccent = (d*point.x)/-(point.z);
    auto yaccent = (d*point.y)/-(point.z);
    return Point2D(xAccent, yaccent);
};

Lines2D doProjection(const Figures3D &f){
    double d = 1;
    Lines2D lines;
    for(auto& fig :f){
        auto faces = fig.faces;
        for (auto face: faces){
            int point1 = face.point_indexes[0];
            int point2 = face.point_indexes[1];

            Vector3D vec1 = fig.points[point1];
            Vector3D vec2 = fig.points[point2];

            auto p1 = doProjectionPoint(vec1, d);
            auto p2 = doProjectionPoint(vec2, d);

            Line2D l(p1, p2, fig.color);
            lines.push_back(l);
        }
    }
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
    auto type = configuration["General"]["type"].as_string_or_die();
    auto size = configuration["General"]["size"].as_int_or_die();

    //achtergrond kleur maken
    auto backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    std::map<int, double> bgcolormap;
    int bgcm = 0;
    for (auto c: backgroundcolor) {
        bgcolormap[bgcm] = c;
        bgcm++;
    }
    Color bgcolor(bgcolormap[0], bgcolormap[1], bgcolormap[2]);

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

        double angle_rad = (angle * M_PI)/180;
        double startingangle_rad = (starting_angle * M_PI)/180;
        double tempangle;

        int iter = 0;
        std::string full_string = initiator;
        while (iter < iterations) { //Terug veranderen naar iterations!!
            std::string temp;
            for (auto &c: full_string) {
                if (c != '+' && c != '-' && c != '(' && c != ')') {
                    auto replacement = lsystem.get_replacement(c);
                    temp += replacement;
                } else{
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
                }
                else{
                    x += std::cos(startingangle_rad);
                    y += std::sin(startingangle_rad);
                }
            } else if (c == '+') {
                startingangle_rad += angle_rad;
            } else if (c == '-') {
                startingangle_rad -= angle_rad;
            }
            else if (c == '('){
                pointend.x = x; pointend.y = y;
                angle_stack.push_back(startingangle_rad);
                bracket_stack.push_back(pointend);
            }
            else if (c == ')'){
                Point2D last_point = bracket_stack.back();
                startingangle_rad = angle_stack.back();
                x = last_point.x; y = last_point.y;
                bracket_stack.pop_back();
                angle_stack.pop_back();
            }
        }
        auto image = Draw2DLines(lijntjes, size, bgcolor, color2);
        return image;
    }

    else if (type == "Wireframe"){
        Figures3D figures;
        auto nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        auto eye = configuration["General"]["eye"].as_double_tuple_or_die();
        for (int i = 0; i < nrFigures; i++){
            std::string figure_name = "figure" + std::to_string(i);
            Figure fig;
            auto secondary_type = configuration[figure_name]["type"].as_string_or_die();
            auto rotateX = configuration[figure_name]["rotateX"].as_double_or_die();
            auto rotateY = configuration[figure_name]["rotateY"].as_double_or_die();
            auto rotateZ = configuration[figure_name]["rotateZ"].as_double_or_die();
            auto scale = configuration[figure_name]["scale"].as_double_or_die();

            auto center = configuration[figure_name]["center"].as_double_tuple_or_die();
            auto drawing_color = configuration[figure_name]["color"].as_double_tuple_or_die();
            std::map<int, double> colormap; int cm = 0;
            for (auto c: drawing_color){
                colormap[cm] = c;
                cm++;
            }
            Color drawingcolor = Color(colormap[0], colormap[1], colormap[2]);

            auto nrpoints = configuration[figure_name]["nrPoints"].as_int_or_die();
            auto nrlines = configuration[figure_name]["nrLines"].as_int_or_die();

            for (int j = 0; j < nrpoints; j++){
                std::string point_name = "point"+ std::to_string(i);
                Vector3D point;
                std::map<int, double> pointmap;
                int pp = 0;
                auto point_tuple = configuration[figure_name][point_name].as_double_tuple_or_die();
                for (auto p: point_tuple){
                    pointmap[pp] = p;
                    pp++;
                }
                point.x = pointmap[0]; point.y = pointmap[1]; point.z = pointmap[2];
                fig.points.push_back(point);
            }
            for (int j = 0; j < nrlines; j++){
                std::string line_name = "line"+ std::to_string(j);
                Face face;
                auto linetuple = configuration[figure_name][line_name].as_int_tuple_or_die();
                for (auto l: linetuple){
                    face.point_indexes.push_back(l);
                }
                fig.faces.push_back(face);
            }
            fig.color = drawingcolor;
            figures.push_back(fig);
        }

    }
    return img::EasyImage();
}


int main(int argc, char const *argv[]) {
    std::ifstream inputfile;
    inputfile.open("line_drawings002.ini");
    ini::Configuration config;
    inputfile >> config;
    inputfile.close();
    generate_image(config);
}
