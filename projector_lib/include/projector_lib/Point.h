#ifndef _POINT_H_
#define _POINT_H_
class Point {
public:
	int x, y;
	Point(int i_x, int i_y) : x(i_x), y(i_y){};
	Point() : x(-999), y(-999) {};
};
#endif