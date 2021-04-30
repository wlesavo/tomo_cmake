#ifndef _POINT_H_
#define _POINT_H_
class Point {
public:
	double x, y, z;
	Point(double i_x, double i_y, double i_z = 0) : x(i_x), y(i_y), z(i_z){};
	Point() : x(-999), y(-999), z(0) {};
};
#endif