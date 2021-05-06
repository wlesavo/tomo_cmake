#ifndef _LINE_HPP_
#define _LINE_HPP_
#include <vector>
#include <projector_lib/Point.h>

class Line {
public:
	bool transpose = false, reverse_x = false, reverse_y = false, reverse_z = false;
	bool isVertical = false, isHorizontal = false;
	double angle = 0.0, k = 0.0, b = 0.0;
	double kx, ky, kz, x0, y0, z0;
	double reverse_k, reverse_ab;

	Line(double x1, double y1, double z1, double x2, double y2, double z2);
	Line(double x1, double y1, double x2, double y2);
	Line(double x, double y, double angle);
	Line() {};
	Line(const Line& line1);

	double value(double x) const;
	double coor(double y) const;
	double getLambda(double position, int axis = 0) const;
	Point getPoint(double lambda) const;

};

#endif