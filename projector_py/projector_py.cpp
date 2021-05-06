#include <projector_lib/Projector.hpp>
#include "projector_py.h"
#include <projector_lib/Geometries.h>
#include <projector_lib/SumAlgo.h>

// geometry classes and destructors

PROJECTOR_API struct GeometryFanBeamPtr* geometry_fan_beam_create(double* i_angles, size_t angles_size, int i_detectorCount, double i_detectorSize, double i_dSourceToObj
	, double i_dObjToDetector, double imgCenterX, double imgCenterY) {
	std::vector<double> angles(i_angles, i_angles + angles_size);
	return reinterpret_cast<struct GeometryFanBeamPtr*> (new std::shared_ptr<GeometryFanBeam>(new GeometryFanBeam(angles, i_detectorCount, i_detectorSize, i_dSourceToObj
		, i_dObjToDetector, imgCenterX, imgCenterY)));
}

PROJECTOR_API void geometry_fan_beam_free(struct GeometryFanBeamPtr* geometry_fan_beam){
    delete reinterpret_cast<std::shared_ptr<GeometryFanBeam>*>(geometry_fan_beam);
}

PROJECTOR_API struct GeometryParallelPtr* geometry_parallel_create(double* i_angles, size_t angles_size, int i_detectorCount, double i_detectorSize, double imgCenterX, double imgCenterY) {
	std::vector<double> angles(i_angles, i_angles + angles_size);
	return reinterpret_cast<struct GeometryParallelPtr*> (new std::shared_ptr<GeometryParallel>(new GeometryParallel(angles, i_detectorCount, i_detectorSize
		, imgCenterX, imgCenterY)));
}

PROJECTOR_API void geometry_parallel_free(struct GeometryParallelPtr* geometry_parallel) {
	delete reinterpret_cast<std::shared_ptr<GeometryFanBeam>*>(geometry_parallel);
}

// projector class and destructor

PROJECTOR_API struct Projector* projector_create(float* inputImg, int imgSize_x, int imgSize_y
	, void* geometry, int sumAlgorithm) {
	size_t size = static_cast<size_t>(imgSize_x) * imgSize_y;
	std::unique_ptr<float[]> img(new float[size]);
	std::copy(inputImg, inputImg + size, img.get());
	return new Projector(std::move(img), *reinterpret_cast<std::shared_ptr<Geometry>*>(geometry), reinterpret_cast<SumAlgorithm&>(sumAlgorithm), imgSize_x, imgSize_y);
}

PROJECTOR_API void projector_free(struct Projector* projector) {
	delete reinterpret_cast<Projector*>(projector);
}

// projector methods

PROJECTOR_API void projector_get_full_projection(float** out, size_t* out_size_x, size_t* out_size_y, struct Projector* projector) {
	*out = projector->getFullProjection().release();
	*out_size_x = static_cast<size_t>(projector->geometry->nDetectors);
	*out_size_y = static_cast<size_t>(projector->geometry->angles.size());
}

PROJECTOR_API void projector_get_single_projection(float** out, size_t* out_size, struct Projector* projector, int angle) {
	*out = projector->getSingleProjection(angle).release();
	*out_size = static_cast<size_t>(projector->geometry->nDetectors);
}

PROJECTOR_API void float_free(float* f) {
	delete f;
}