#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#include <vector>
#include <cmath>
class Geometry {
public:
	const double dSourceToObj, dObjToDetector, detectorSize, imgCenterX, imgCenterY, imgCenterZ;
	const int detectorCount, detectorCount_x, detectorCount_y;
	const std::vector<double> angles;
	double part;
	std::vector<double> distanceToDetectorPixelCenter,	angleToDetectorPixelCenter;
	std::vector<double> distanceToDetectorPixelLeft,	angleToDetectorPixelLeft;
	std::vector<double> distanceToDetectorPixelRight,	angleToDetectorPixelRight;
	std::vector<double> distanceToDetectorPixel,		angleToDetectorPixel;

	Geometry(const std::vector<double>& i_angles, int i_detectorCount_x, int i_detectorCount_y, double i_detectorSize, double i_dSourceToObj
		, double i_dObjToDetector, double imgCenterX, double imgCenterY, double imgCenterZ)
		: dSourceToObj(i_dSourceToObj), dObjToDetector(i_dObjToDetector), detectorSize(i_detectorSize)
		, detectorCount(i_detectorCount_x * i_detectorCount_y), detectorCount_x(i_detectorCount_x), detectorCount_y(i_detectorCount_y)
		, angles(i_angles), imgCenterX(imgCenterX), imgCenterY(imgCenterY), imgCenterZ(imgCenterZ)
	{};
	
	virtual Line v_GetNextLineCenter(int angle, int detector_i) = 0;
	virtual std::pair<Line,Line> v_GetNextLinePair(int angle, int detector_i, bool parallel = true) = 0;
	virtual ~Geometry() {};
	
};
#endif