#ifndef _PROJECTOR_HPP_
#define _PROJECTOR_HPP_
#include <projector_lib/Line.h>
#include <projector_lib/MyImg.h>
#include <projector_lib/Geometry.h>
#include <projector_lib/Point.h>
#include <projector_lib/SumAlgorithm.h>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

class Projector {
public:
	// main objects
	const MyImg inputImg;
	SumAlgorithm sumAlgorithm;

	const int imgSize_x, imgSize_y, imgSize_z;
	std::shared_ptr<Geometry> geometry;
	std::unique_ptr<float[]> fullProjection = nullptr;

	// debug vars
	bool testFlag = false;

	// constructor
	// todo: sizes to geometry
	Projector(std::unique_ptr<float[]> inputImg, std::shared_ptr<Geometry> geometry, SumAlgorithm sumAlgorithm, int imgSize_x, int imgSize_y, int imgSize_z = 1);
	
	// main methods
	std::unique_ptr<float[]> getSingleProjection(int angle_i) const;
	void buildFullProjection();
	std::unique_ptr<float[]> getFullProjection();
	std::unique_ptr<unsigned char[]> getFullProjectionImage();

	// common methods
	std::pair<Point, Point> getIntersectionPoints(const Line& line) const;
	Line constructLine(const Line& line) const;
	Line constructLine(const Line& line, bool transpose, bool reverse_x) const;
	Line constructLine3D(const Line& line) const;
	std::pair<Line, Line> sortLines(const Line& line1, const Line& line2) const;

	float manyPixelArea(int i_min, int i_max, int j, bool isUpper, const Line& line) const;
	double singlePixelArea(int i, int j, const Line& line) const;
	float sumNeibs(double i_min, double i_max, double j, bool transpose, bool reverse_x, int slice = 0) const;

	// sum algorithms
	float sumLine(const Line& line, int slice = 0) const;
	float sumLine3D(const Line& line) const;
	float sumLinear(const Line& line) const;
	float sumArea(const Line& line1, const Line& line2) const;
	float sumAreaExact(const Line& line_1, const Line& line_2) const;
	
	// debug methods
	float getLineProjectionTest(int angle, int detector);

};

#endif