#include <projector_lib/Line.h>
#include <cmath>
#include <cassert>

Line::Line(double x1, double y1, double x2, double y2){
	angle = std::atan2(y2 - y1, x2 - x1);
	if (std::abs(std::cos(angle)) < 1e-10) {
		vert = true;
		b = x1;
	}
	else if (std::abs(std::sin(angle)) < 1e-10) {
		hor = true;
		k = 0;
		b = y1;
	}
	else {
		k = (y2 - y1) / (x2 - x1);// std::tan(angle);
		b = y2 - k * x2;
		reverse_k = 1.0 / k;
		reverse_ab = 1 / std::pow(1 + k * k, 0.5);
	}
}

Line::Line(double x, double y, double i_angle){
	angle = i_angle;
	if (std::abs(std::cos(angle)) < 1e-10) {
		vert = true;
		b = x;
	}
	else if (std::abs(std::sin(angle)) < 1e-10) {
		hor = true;
		k = 0;
		b = y;
	}
	else {
		k = std::tan(angle);
		b = y - k * x;
		reverse_k = 1.0 / k;
		reverse_ab = 1 / std::pow(1 + k * k, 0.5);
	}
}

Line::Line(double x1, double y1, double z1, double x2, double y2, double z2) {
	double d, dx, dy, dz;
	dx = x2 - x1;
	dy = y2 - y1;
	dz = z2 - z1;
	d = std::pow(dx * dx + dy * dy + dz * dz, 0.5);
	assert(d > 0);
	kx = dx / d;
	ky = dy / d;
	kz = dz / d;
	x0 = x1;
	y0 = y1;
	z0 = z1;
}

Line::Line(const Line& line) : transpose(line.transpose), reverse_x(line.reverse_x), reverse_y(line.reverse_y), reverse_z(line.reverse_z)
	, vert(line.vert), hor(line.hor), angle(line.angle), k(line.k), b(line.b), reverse_k(line.reverse_k), reverse_ab(line.reverse_ab),
	kx(line.kx), ky(line.ky), kz(line.kz), x0(line.x0), y0(line.y0), z0(line.z0)
{}

double Line::value(double x) const {
	return x * k + b;
}

double Line::coor(double y) const {
	//double rev = 1 / k ;
	return ((y - b) / k);
};

double Line::getLambda(double pos, int axis) const {
	if (axis == 0) {
		return (pos - x0) / kx;
	}
	if (axis == 1) {
		return (pos - y0) / ky;
	}
	if (axis == 2) {
		return (pos - z0) / kz;
	}
}

Point Line::getPoint(double lambda) const {
	return Point(x0 + lambda * kx, y0 + lambda * ky, z0 + lambda * kz);
}