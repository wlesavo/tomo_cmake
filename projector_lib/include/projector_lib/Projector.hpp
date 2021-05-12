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
	std::unique_ptr<float[]> forwardProjection = nullptr;
	std::unique_ptr<float[]> backProjection = nullptr;

	// debug vars
	bool testFlag = false;

	// constructor
	// todo: sizes to geometry
	Projector(std::unique_ptr<float[]> inputImg, std::shared_ptr<Geometry> geometry, SumAlgorithm sumAlgorithm, int imgSize_x, int imgSize_y, int imgSize_z = 1);
	
	// main methods
	std::unique_ptr<float[]> getForwardProjection();
	std::unique_ptr<unsigned char[]> getForwardProjectionImage();
	
	// weighted methods
	void buildBackProjection();
	void buildForwardProjection();
	void getWeights(int angleIndex, int detectorIndex, int* coorDst, float* weightsDst, int* sizeDst);
	void projectBack(float value, int* coors, float* weights, int size);
	void projectForward(int angleIndex, int detectorIndex, int* coors, float* weights, int size);
	
	// common methods
	std::pair<Point, Point> getIntersectionPoints(const Line& line) const;
	Line constructLine(const Line& line) const;
	Line constructLine(const Line& line, bool transpose, bool reverse_x) const;
	Line constructLine3D(const Line& line) const;
	std::pair<Line, Line> sortLines(const Line& line1, const Line& line2) const;

	float manyPixelArea(int i_min, int i_max, int j, bool isUpper, const Line& line) const;
	double singlePixelArea(int i, int j, const Line& line) const;

	// sum algorithms
	float sumAreaExact(const Line& line_1, const Line& line_2) const;

	// weight algorithms
	void weightNeibsArea(int i_min, int i_max, int j, const Line& line, const Line& line2,
		int* coorDst, float* weightsDst, int* sizeDst) const;
	void weightNeibsLine(double j_min, double j_max, double i, 
		int* coorDst, float* weightsDst, int* sizeDst, bool transpose, bool reverse_x, int slice = 0, float koeff = 1) const;
	
	void getWeightsLine(const Line& line, int* coorDst, float* weightsDst, int* size, int slice = 0, bool isBinary = false) const;
	void getWeightsLine3D(const Line& line, int* coorDst, float* weightsDst, int* sizeDst) const;
	void getWeightsArea(const Line& line_1, const Line& line_2, int* coorDst, float* weightsDst, int* sizeDst) const;
	


};

#endif