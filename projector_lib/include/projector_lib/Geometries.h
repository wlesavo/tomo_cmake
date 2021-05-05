#ifndef _GEOMETRIES_H_
#define _GEOMETRIES_H_
#include <cmath>
#include <vector>
#include <corecrt_math_defines.h>

class GeometryFanBeam : public Geometry {
public:

	GeometryFanBeam(const std::vector<double>& i_angles, int i_detectorCount, double i_detectorSize, double i_dSourceToObj, double i_dObjToDetector, double imgCenterX, double imgCenterY) :
		Geometry(i_angles, i_detectorCount, 1, i_detectorSize, i_dSourceToObj, i_dObjToDetector, imgCenterX, imgCenterY, 0)
	{
		// precalculating distances and angles for pixels of detector matrix
		// in experiment pixels are rotating around origin with new distance and initial angle
		double delta_angle = std::atan2(0.5 - 0.0, i_dSourceToObj + i_dObjToDetector);
		double d1 = (i_dSourceToObj - imgCenterX);
		double d2 = (i_dSourceToObj + imgCenterX);
		part = (1 - (d2 * d2 - d1 * d1) * std::sin(delta_angle) / (imgCenterX * 2)) / 2;

		double center = (detectorCount) * 0.5 * detectorSize;
		for (int i = 0; i < detectorCount; ++i) {
			double det_y = center - detectorSize * (i + 0.5);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixelCenter.push_back(dist);
			angleToDetectorPixelCenter.push_back(angle);
		}
		for (int i = 0; i < detectorCount + 1; ++i) {
			double det_y = center - detectorSize * (i);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixel.push_back(dist);
			angleToDetectorPixel.push_back(angle);
		}

		for (int i = 0; i < detectorCount; ++i) {
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

	virtual Line v_GetNextLineCenter(int i_angle, int detector_i) {
		double angle = M_PI - angles[i_angle];
		double x_source = - dSourceToObj * std::sin(angle) + imgCenterX;
		double y_source = dSourceToObj * std::cos(angle) + imgCenterY;
		double x_det = - distanceToDetectorPixelCenter[detector_i] * std::sin(M_PI + angle + angleToDetectorPixelCenter[detector_i]) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detector_i] * std::cos(M_PI + angle + angleToDetectorPixelCenter[detector_i]) + imgCenterY;
		return Line(x_source, y_source, x_det, y_det);
	};

	virtual std::pair<Line, Line> v_GetNextLinePair(int i_angle, int detector_i, bool parallel = true) {
		double angle = angles[i_angle];
		if (parallel) {
			Line line = v_GetNextLineCenter(i_angle, detector_i);
			double x_det1 = distanceToDetectorPixelLeft[detector_i] * std::sin(angle - angleToDetectorPixelLeft[detector_i]) + imgCenterX;
			double y_det1 = distanceToDetectorPixelLeft[detector_i] * std::cos(angle - angleToDetectorPixelLeft[detector_i]) + imgCenterY;
			double x_det2 = distanceToDetectorPixelRight[detector_i] * std::sin(angle - angleToDetectorPixelRight[detector_i]) + imgCenterX;
			double y_det2 = distanceToDetectorPixelRight[detector_i] * std::cos(angle - angleToDetectorPixelRight[detector_i]) + imgCenterY;
			return std::pair<Line, Line>(Line(x_det1, y_det1, line.angle), Line(x_det2, y_det2, line.angle));
		}
		if (!parallel) {
			double x_source = -dSourceToObj * std::sin(angle) + imgCenterX;
			double y_source = -dSourceToObj * std::cos(angle) + imgCenterY;
			double x_det1 = distanceToDetectorPixel[detector_i] * std::sin(angle - angleToDetectorPixel[detector_i]) + imgCenterX;
			double y_det1 = distanceToDetectorPixel[detector_i] * std::cos(angle - angleToDetectorPixel[detector_i]) + imgCenterY;
			double x_det2 = distanceToDetectorPixel[detector_i+1] * std::sin(angle - angleToDetectorPixel[detector_i+1]) + imgCenterX;
			double y_det2 = distanceToDetectorPixel[detector_i+1] * std::cos(angle - angleToDetectorPixel[detector_i+1]) + imgCenterY;
			return std::pair<Line, Line>(Line(x_source, y_source, x_det1, y_det1), Line(x_source, y_source, x_det2, y_det2));
		}
	};

};

class GeometryParallel : public Geometry {
public:
	GeometryParallel(const std::vector<double>& i_angles, int i_detectorCount, double i_detectorSize, double imgCenterX, double imgCenterY) :
		Geometry(i_angles, i_detectorCount, 1, i_detectorSize, 0, 0, imgCenterX, imgCenterY, 0)
	{
		// precalculating distances and angles for pixels of detector matrix
		// in experiment pixels are rotating around origin with new distance and initial angle

		double center = (detectorCount) * 0.5 * detectorSize;
		for (int i = 0; i < detectorCount; ++i) {
			double det_y = detectorSize * (i + 0.5 - center);
			distanceToDetectorPixelCenter.push_back(det_y);
		}

		for (int i = 0; i < detectorCount + 1; ++i) {
			double det_y = detectorSize * (i - center);
			distanceToDetectorPixel.push_back(det_y);
		}
	}

	virtual Line v_GetNextLineCenter(int i_angle, int detector_i) {
		double angle = angles[i_angle];
		double x_det = distanceToDetectorPixelCenter[detector_i] * std::cos(-angle) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detector_i] * std::sin(-angle) + imgCenterY;
		return Line(x_det, y_det, M_PI_2 - angle);
	};

	virtual std::pair<Line, Line> v_GetNextLinePair(int i_angle, int detector_i, bool parallel = true) {
		double angle = angles[i_angle];
		double x_det1 = distanceToDetectorPixel[detector_i] * std::cos(-angle) + imgCenterX;
		double y_det1 = distanceToDetectorPixel[detector_i] * std::sin(-angle) + imgCenterY;
		double x_det2 = distanceToDetectorPixel[detector_i+1] * std::cos(-angle) + imgCenterX;
		double y_det2 = distanceToDetectorPixel[detector_i+1] * std::sin(-angle) + imgCenterY;
		return std::pair<Line, Line>(Line(x_det1, y_det1, M_PI_2 - angle), Line(x_det2, y_det2, M_PI_2 - angle));
	};

};

class GeometryFanBeam3D : public Geometry {
public:

	GeometryFanBeam3D(const std::vector<double>& i_angles, int i_detectorCount_x, int i_detectorCount_y, double i_detectorSize, double i_dSourceToObj, double i_dObjToDetector, double imgCenterX, double imgCenterY, double imgCenterZ) :
		Geometry(i_angles, i_detectorCount_x, i_detectorCount_y, i_detectorSize, i_dSourceToObj, i_dObjToDetector, imgCenterX, imgCenterY, imgCenterZ)
	{
		// precalculating distances and angles for pixels of detector matrix
		// in experiment pixels are rotating around origin with new distance and initial angle
		
		double center = (detectorCount_x) * 0.5 * detectorSize;
		for (int i = 0; i < detectorCount_x; ++i) {
			double det_y = center - detectorSize * (i + 0.5);
			double dist = std::pow(dObjToDetector * dObjToDetector + det_y * det_y, 0.5);
			double angle = std::atan2(-det_y, dObjToDetector);
			distanceToDetectorPixelCenter.push_back(dist);
			angleToDetectorPixelCenter.push_back(angle);
		}

	}

	virtual Line v_GetNextLineCenter(int i_angle, int detector) {
		int detector_i = detector % detectorCount_x;
		int detector_j = detector / detectorCount_x;
		double angle = M_PI - angles[i_angle];
		double x_source = -dSourceToObj * std::sin(angle) + imgCenterX;
		double y_source = dSourceToObj * std::cos(angle) + imgCenterY;
		double z_source = imgCenterZ;
		double x_det = -distanceToDetectorPixelCenter[detector_i] * std::sin(M_PI + angle + angleToDetectorPixelCenter[detector_i]) + imgCenterX;
		double y_det = distanceToDetectorPixelCenter[detector_i] * std::cos(M_PI + angle + angleToDetectorPixelCenter[detector_i]) + imgCenterY;
		double z_det = detectorSize * (detectorCount_y * 0.5 - (detector_j + 1)) + imgCenterZ;
		return Line(x_source, y_source, z_source, x_det, y_det, z_det);
	};

	virtual std::pair<Line, Line> v_GetNextLinePair(int i_angle, int detector_i, bool parallel = true) {
		return std::pair<Line, Line>(Line(), Line());		
	};
};


#endif