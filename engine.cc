#include "easy_image.h"
#include "ini_configuration.h"
#include "LSystems.h"
#include "l_parser/l_parser.h"

#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <set>
#include <cmath>

using Lines2D = std::list<Line2D>;

img::EasyImage Draw2DLines(Lines2D &lines, const int size, Color bgc, Color drawing_color) {
    auto backgroundcolor = img::Color(bgc.red*255, bgc.green*255, bgc.blue*255);
    auto drawingcolor = img::Color(drawing_color.red*255, drawing_color.green*255, drawing_color.blue*255);

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
    std::list<Point2D> all_points;
    for (auto& l: lines){
        auto points = l.get_coordinates();
        for (auto p: points){
            all_points.push_back(p);
        }
    }
    for (auto p: all_points){
        image(p.x, p.y) = drawingcolor;
    }
    return image;
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
    return img::EasyImage();
}


int main(int argc, char const *argv[]) {
    ini::Configuration input;
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
