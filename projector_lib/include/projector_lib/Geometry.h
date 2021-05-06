#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#include <vector>
#include <cmath>
class Geometry {
public:
	const double dSourceToObj, dObjToDetector, detectorSize, imgCenterX, imgCenterY, imgCenterZ;
	const int nDetectors, nDetectorsX, nDetectorsY;
	const std::vector<double> angles;
	std::vector<double> distanceToDetectorPixelCenter, angleToDetectorPixelCenter;
	std::vector<double> distanceToDetectorPixel,	   angleToDetectorPixel;

	Geometry(const std::vector<double>& angles, int nDetectorsX, int nDetectorsY, double detectorSize, double dSourceToObj
		, double dObjToDetector, double imgCenterX, double imgCenterY, double imgCenterZ)
		: dSourceToObj(dSourceToObj), dObjToDetector(dObjToDetector), detectorSize(detectorSize)
		, nDetectors(nDetectorsX * nDetectorsY), nDetectorsX(nDetectorsX), nDetectorsY(nDetectorsY)
		, angles(angles), imgCenterX(imgCenterX), imgCenterY(imgCenterY), imgCenterZ(imgCenterZ)
	{};
	
	virtual Line getLine(int angleIndex, int detectorIndex) = 0;
	virtual std::pair<Line,Line> getLinePair(int angleIndex, int detectorIndex, bool isParallel = true) = 0;
	virtual ~Geometry() {};
	
};
#endif