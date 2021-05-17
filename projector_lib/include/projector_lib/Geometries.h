#ifndef _GEOMETRIES_H_
#define _GEOMETRIES_H_
#include <cmath>
#include <vector>
#include <corecrt_math_defines.h>

class GeometryFanBeam : public Geometry {
public:

	std::vector<double> distanceToDetectorPixelLeft, angleToDetectorPixelLeft;
	std::vector<double> distanceToDetectorPixelRight, angleToDetectorPixelRight;

	GeometryFanBeam(const std::vector<double>& angles, int nDetectors, double detectorSize, double dSourceToObj, double dObjToDetector, double imgCenterX, double imgCenterY) :
		Geometry(angles, nDetectors, 1, detectorSize, dSourceToObj, dObjToDetector, imgCenterX, imgCenterY, 0)
	{
		// Precalculating distances and angles for pixels of detector matrix.
		// In the experiment pixels are rotating around origin with new distance and initial angle.
		// Vectors Left and Right are used for approximation of parallel beam with the same area
		// as original fan beam.

		double delta_angle = std::atan2(0.5 - 0.0, dSourceToObj + dObjToDetector);
		double d1 = (dSourceToObj - imgCenterX);
		double d2 = (dSourceToObj + imgCenterX);
		double part = (1 - (d2 * d2 - d1 * d1) * std::sin(delta_angle) / (imgCenterX * 2)) / 2;

		double center = (nDetectors) * 0.5 * detectorSize;
		
		for (int i = 0; i < nDetectors; ++i) {
			double det_y = center - detectorSize * (i + 0.5);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixelCenter.push_back(dist);
			angleToDetectorPixelCenter.push_back(angle);
		}
		
		for (int i = 0; i < nDetectors + 1; ++i) {
			double det_y = center - detectorSize * (i);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixel.push_back(dist);
			angleToDetectorPixel.push_back(angle);
		}

		for (int i = 0; i < nDetectors; ++i) {
			double det_y = center - detectorSize * (i + part);
			//double det_y = center - detectorSize * (i);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(det_y, dObjToDetector);
			distanceToDetectorPixelLeft.push_back(dist);
			angleToDetectorPixelLeft.push_back(angle);
			det_y = center - detectorSize * (i + 1 - part);
			//det_y = center - detectorSize * (i + 1);
			dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			angle = std::atan2(det_y, dObjToDetector);
			distanceToDetectorPixelRight.push_back(dist);
			angleToDetectorPixelRight.push_back(angle);
		}

	}

	virtual Line getLine(int angleIndex, int detectorIndex) {

		double angle = M_PI - angles[angleIndex];
		double x_source = - dSourceToObj * std::sin(angle) + imgCenterX;
		double y_source =   dSourceToObj * std::cos(angle) + imgCenterY;
		double x_det = - distanceToDetectorPixelCenter[detectorIndex] * std::sin(M_PI + angle + angleToDetectorPixelCenter[detectorIndex]) + imgCenterX;
		double y_det =	 distanceToDetectorPixelCenter[detectorIndex] * std::cos(M_PI + angle + angleToDetectorPixelCenter[detectorIndex]) + imgCenterY;
		return Line(x_source, y_source, x_det, y_det);
	};

	virtual std::pair<Line, Line> getLinePair(int angleIndex, int detectorIndex, bool isParallel = true) {
		double angle = angles[angleIndex];
		if (!isParallel) {
			double x_source = - dSourceToObj * std::sin(angle) + imgCenterX;
			double y_source = - dSourceToObj * std::cos(angle) + imgCenterY;
			double x_det1 = distanceToDetectorPixel[detectorIndex] * std::sin(angle - angleToDetectorPixel[detectorIndex]) + imgCenterX;
			double y_det1 = distanceToDetectorPixel[detectorIndex] * std::cos(angle - angleToDetectorPixel[detectorIndex]) + imgCenterY;
			double x_det2 = distanceToDetectorPixel[detectorIndex + 1] * std::sin(angle - angleToDetectorPixel[detectorIndex + 1]) + imgCenterX;
			double y_det2 = distanceToDetectorPixel[detectorIndex + 1] * std::cos(angle - angleToDetectorPixel[detectorIndex + 1]) + imgCenterY;
			return std::pair<Line, Line>(Line(x_source, y_source, x_det1, y_det1), Line(x_source, y_source, x_det2, y_det2));
		}
		else if (isParallel) {
			Line line = getLine(angleIndex, detectorIndex);
			double x_det1 = distanceToDetectorPixelLeft[detectorIndex] * std::sin(angle - angleToDetectorPixelLeft[detectorIndex]) + imgCenterX;
			double y_det1 = distanceToDetectorPixelLeft[detectorIndex] * std::cos(angle - angleToDetectorPixelLeft[detectorIndex]) + imgCenterY;
			double x_det2 = distanceToDetectorPixelRight[detectorIndex] * std::sin(angle - angleToDetectorPixelRight[detectorIndex]) + imgCenterX;
			double y_det2 = distanceToDetectorPixelRight[detectorIndex] * std::cos(angle - angleToDetectorPixelRight[detectorIndex]) + imgCenterY;
			return std::pair<Line, Line>(Line(x_det1, y_det1, line.angle), Line(x_det2, y_det2, line.angle));
		}
	};

};

class GeometryParallel : public Geometry {
public:
	GeometryParallel(const std::vector<double>& angles, int nDetectors, double detectorSize, double imgCenterX, double imgCenterY) :
		Geometry(angles, nDetectors, 1, detectorSize, 0, 0, imgCenterX, imgCenterY, 0)
	{
		// Precalculating distances and angles for pixels of detector matrix.
		// In the experiment pixels are rotating around origin with new distance and initial angle.

		double center = (nDetectors) * 0.5 * detectorSize;
		for (int i = 0; i < nDetectors; ++i) {
			double det_y = detectorSize * (i + 0.5 - center);
			distanceToDetectorPixelCenter.push_back(det_y);
		}

		for (int i = 0; i < nDetectors + 1; ++i) {
			double det_y = detectorSize * (i - center);
			distanceToDetectorPixel.push_back(det_y);
		}
	}

	virtual Line getLine(int angleIndex, int detectorIndex) {
		double angle = angles[angleIndex];
		double x_det = distanceToDetectorPixelCenter[detectorIndex] * std::cos(-angle) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detectorIndex] * std::sin(-angle) + imgCenterY;
		return Line(x_det, y_det, M_PI_2 - angle);
	};

	virtual std::pair<Line, Line> getLinePair(int angleIndex, int detectorIndex, bool isParallel = true) {
		double angle = angles[angleIndex];
		double x_det1 = distanceToDetectorPixel[detectorIndex] * std::cos(-angle) + imgCenterX;
		double y_det1 = distanceToDetectorPixel[detectorIndex] * std::sin(-angle) + imgCenterY;
		double x_det2 = distanceToDetectorPixel[detectorIndex + 1] * std::cos(-angle) + imgCenterX;
		double y_det2 = distanceToDetectorPixel[detectorIndex + 1] * std::sin(-angle) + imgCenterY;
		return std::pair<Line, Line>(Line(x_det1, y_det1, M_PI_2 - angle), Line(x_det2, y_det2, M_PI_2 - angle));
	};

};

class GeometryFanBeam3D : public Geometry {
public:

	GeometryFanBeam3D(const std::vector<double>& angles, int nDetectorsX, int nDetectorsY, double detectorSize, double dSourceToObj, double dObjToDetector, double imgCenterX, double imgCenterY, double imgCenterZ) :
		Geometry(angles, nDetectorsX, nDetectorsY, detectorSize, dSourceToObj, dObjToDetector, imgCenterX, imgCenterY, imgCenterZ)
	{
		// Precalculating distances and angles for pixels of detector matrix.
		// In the experiment pixels are rotating around origin with new distance and initial angle.

		double center = (nDetectorsX) * 0.5 * detectorSize;
		for (int i = 0; i < nDetectorsX; ++i) {
			double det_y = center - detectorSize * (i + 0.5);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixelCenter.push_back(dist);
			angleToDetectorPixelCenter.push_back(angle);
		}

	}

	virtual Line getLine(int angleIndex, int detectorIndex) {
		int detector_i = detectorIndex % nDetectorsX;
		int detector_j = detectorIndex / nDetectorsX;
		double angle = M_PI - angles[angleIndex];
		double x_source = - dSourceToObj * std::sin(angle) + imgCenterX;
		double y_source =   dSourceToObj * std::cos(angle) + imgCenterY;
		double z_source = imgCenterZ;
		double x_det = - distanceToDetectorPixelCenter[detector_i] * std::sin(M_PI + angle + angleToDetectorPixelCenter[detector_i]) + imgCenterX;
		double y_det =   distanceToDetectorPixelCenter[detector_i] * std::cos(M_PI + angle + angleToDetectorPixelCenter[detector_i]) + imgCenterY;
		double z_det = detectorSize * (nDetectorsY * 0.5 - (detector_j + 0.5)) + imgCenterZ;
		return Line(x_source, y_source, z_source, x_det, y_det, z_det);
	};

	virtual std::pair<Line, Line> getLinePair(int angleIndex, int detectorIndex, bool isParallel = true) {
		return std::pair<Line, Line>(Line(), Line());		
	};
};

class GeometryParallel3D : public Geometry {
public:

	GeometryParallel3D(const std::vector<double>& angles, int nDetectorsX, int nDetectorsY, double detectorSize, double imgCenterX, double imgCenterY, double imgCenterZ) :
		Geometry(angles, nDetectorsX, nDetectorsY, detectorSize, 0.0, 0.0, imgCenterX, imgCenterY, imgCenterZ)
	{
		// Precalculating distances and angles for pixels of detector matrix.
		// In the experiment pixels are rotating around origin with new distance and initial angle.

		double center = (nDetectorsX) * 0.5 * detectorSize;
		for (int i = 0; i < nDetectorsX; ++i) {
			double det_y = detectorSize * ((i + 0.5) - center);
			distanceToDetectorPixelCenter.push_back(det_y);
		}
	}

	virtual Line getLine(int angleIndex, int detectorIndex) {
		double angle = - angles[angleIndex];
		double x_det = distanceToDetectorPixelCenter[detectorIndex] * std::cos(-angle) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detectorIndex] * std::sin(-angle) + imgCenterY;
		return Line(x_det, y_det, M_PI_2 - angle);
	};

	virtual std::pair<Line, Line> getLinePair(int angleIndex, int detectorIndex, bool isParallel = true) {
		return std::pair<Line, Line>(Line(), Line());
	};

};

#endif