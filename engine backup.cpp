//
// Created by olivi on 3/17/2023.
//
img::EasyImage Draw2DLines(Lines2D &lines, const int size, Color bgc, Color drawing_color) {
    auto background = img::Color(bgc.red*255, bgc.green*255, bgc.blue*255);
    double xmin = size, xmax = 0, ymin = size, ymax = 0;

    for (auto& l: lines) {
        if (l.p1.x < xmin) xmin = l.p1.x;
        else if (l.p1.x > xmax) xmax = l.p1.x;
        if (l.p1.y < ymin) ymin = l.p1.y;
        else if (l.p1.y > ymax) ymax = l.p1.y;
        if (l.p2.x < xmin) xmin = l.p2.x;
        else if (l.p2.x > xmax) xmax = l.p2.x;
        if (l.p2.y < ymin) ymin = l.p2.y;
        else if (l.p2.y > ymax) ymax = l.p2.y;

    }
    //grootte image berekenen
    auto xrange = xmax - xmin;
    auto yrange = ymax - ymin;
    unsigned int imagex = size * (xrange / std::max(xrange, yrange));
    unsigned int imagey = size * (yrange / std::max(xrange, yrange));

    //schaalfactor berekenen
    auto schaalfactor_d = 0.95 * (imagex / xrange);
    for (auto &l: lines) {
        l.p1.x = l.p1.x * schaalfactor_d;
        l.p1.y = l.p1.y * schaalfactor_d;
        l.p2.x = l.p2.x * schaalfactor_d;
        l.p2.y = l.p2.y * schaalfactor_d;
    }

    //tekening verschuiven
    double DCx = schaalfactor_d * ((xmax + xmin) / 2);
    double DCy = schaalfactor_d * ((ymax + ymin) / 2);
    double dx = imagex / 2 - DCx;
    double dy = imagey / 2 - DCy;

    for (auto &l: lines) {
        l.p1.x = std::round(l.p1.x + dx);
        l.p1.y = std::round(l.p1.y + dy);
        l.p2.x = std::round(l.p2.x + dx);
        l.p2.y = std::round(l.p2.y + dy);
    }
    auto drawingcolor = img::Color(drawing_color.red*255,drawing_color.green*255,drawing_color.blue*255);
    std::list<Point2D> all_points;
    for (auto& l: lines){
        auto line_points = l.get_coordinates();
        for (auto i: line_points){
            all_points.push_back(i);
        }
    }
    img::EasyImage image(imagex, imagey, background);
    for (auto p: all_points) {
        image(p.x, p.y) = drawingcolor ;
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
                bracket_stack.push_back(pointend);
            }
            else if (c == ')'){
                Point2D last_point = bracket_stack.back();
                x = last_point.x; y = last_point.y;
                bracket_stack.pop_back();
            }
        }
        auto image = Draw2DLines(lijntjes, size, bgcolor, color2);
        return image;
    }
    return img::EasyImage();
}