#ifndef _POINT_H_
#define _POINT_H_
class Point {
public:
	double x, y, z;
	Point(double x, double y, double z = 0.0) : x(x), y(y), z(z){};
	Point() : x(-999.0), y(-999.0), z(0.0) {};
};
#endif