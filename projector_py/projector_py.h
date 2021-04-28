#if defined _MSC_VER
# define PROJECTOR_API __declspec(dllexport)
#elif defined(__clang__) || defined(__GNUC__)
# define PROJECTOR_API  \
  __attribute__ ((visibility ("default")))  \
  __attribute__ ((noinline))
#else
# define PROJECTOR_API
#endif

#ifdef __cplusplus
extern "C"
{
#endif
// geometry classes and destructors

PROJECTOR_API struct GeometryFanBeamPtr* geometry_fan_beam_create(double* i_angles, size_t angles_size
	, int i_detectorCount, double i_detectorSize
	, double i_dSourceToObj, double i_dObjToDetector, double imgCenterX, double imgCenterY);
    
PROJECTOR_API void geometry_fan_beam_free(struct GeometryFanBeamPtr* geometry_fan_beam);

PROJECTOR_API struct GeometryParallelPtr* geometry_parallel_create(double* i_angles, size_t angles_size, int i_detectorCount, double i_detectorSize, double imgCenterX, double imgCenterY);

PROJECTOR_API void geometry_parallel_free(struct GeometryParallelPtr* geometry_parallel);


// projector class and destructor

PROJECTOR_API struct Projector* projector_create(float* inputImg, int imgSize_x, int imgSize_y
	, void* geometry, int sumAlgorithm);

PROJECTOR_API void projector_free(struct Projector* projector);

// projector methods

PROJECTOR_API void projector_get_full_projection(float** out, size_t* out_size_x, size_t* out_size_y, struct Projector* projector);

PROJECTOR_API void projector_get_single_projection(float** out, size_t* out_size, struct Projector* projector, int angle);

PROJECTOR_API void float_free(float* f);

#ifdef __cplusplus
}
#endif
