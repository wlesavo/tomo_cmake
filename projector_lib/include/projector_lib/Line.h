#ifndef _LINE_HPP_
#define _LINE_HPP_
#include <vector>

class Line {
public:
	bool transpose = false, reverse_x = false, reverse_y = false;
	bool vert = false, hor = false;
	double angle = 0.0, k = 0.0, b = 0.0;
	double reverse_k, reverse_ab;
	Line(double x1, double y1, double x2, double y2);
	Line(double x, double y, double angle);
	Line() {};
	Line(const Line& line1);
	double value(double x) const;
	double coor(double y) const;
};

#endif