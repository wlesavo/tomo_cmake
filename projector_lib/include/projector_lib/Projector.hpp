#ifndef _PROJECTOR_HPP_
#define _PROJECTOR_HPP_
#include <projector_lib/Line.h>
#include <projector_lib/MyImg.h>
#include <projector_lib/Geometry.h>
#include <projector_lib/Point.h>
#include <projector_lib/SumAlgo.h>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

class Projector {
public:
	// main objects
	const MyImg inputImg;
	SumAlgo sumAlgorithm;

	const int imgSize_x, imgSize_y;
	std::shared_ptr<Geometry> geometry;
	std::unique_ptr<float[]> fullProjection = nullptr;

	// debug vars
	bool testFlag = false;

	// constructor
	Projector(std::unique_ptr<float[]> inputImg, int imgSize_x, int imgSize_y
		, std::shared_ptr<Geometry> geometry, SumAlgo sumAlgorithm);
	
	// main methods

	std::unique_ptr<float[]> getSingleProjection(int angle_i) const;
	void buildFullProjection();
	std::unique_ptr<float[]> getFullProjection();
	std::unique_ptr<unsigned char[]> getFullProjectionImage();

	// common methods
	std::pair<Point, Point> getIntersectionPoints(const Line& line, int pixel_i = -999, int pixel_j = -999) const;
	std::pair<Line, Line> sortLines(const Line& line1, const Line& line2) const;
	float manyPixelArea(int i_min, int i_max, int j, bool upper, const Line& line) const;
	float sumNeibs(double i_min, double i_max, double j, bool transpose, bool reverse_x) const;
	double singlePixelArea(int i, int j, const Line& line) const;
	Line constructLine(const Line& line, bool transpose, bool reverse_x) const;
	Line constructLine(const Line& line) const;

	// sum algorithms
	float sumLine(const Line& line) const;
	float sumLinear(const Line& line) const;
	float sumArea(const Line& line1, const Line& line2) const;
	float sumAreaExact(const Line& line_1, const Line& line_2) const;
	
	// debug methods
	float sumLineTest(const Line& line);
	float getLineProjectionTest(const int& angle, const int& detector);

};

#endif