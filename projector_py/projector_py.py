from cffi import FFI
import re
import numpy as np

def _c_proprocessor(text):
    r_macro = re.compile(r"^\s*#\s*(define|include|pragma|else|endif|ifdef|ifndef|if)\s*(.*)")
    emit = []
    for line in text.splitlines():
        match = r_macro.match(line)
        if match:
            operation = match.group(1)
            if operation in ('pragma', 'include', 'define'):
                pass
            elif operation in ('ifdef', 'if'):
                emit.append(False)
            elif operation == 'ifndef':
                emit.append(True)
            elif operation == 'else':
                emit[-1] = not emit[-1]
            elif operation == 'endif':
                emit.pop()
            else:
                assert False
            continue
        if all(emit):
            yield line
    assert not emit


ffi = FFI()
with open("projector_py.h") as f:
    s = "\n".join(_c_proprocessor(f.read()))
    s = s.replace("PROJECTOR_API", "")
    ffi.cdef(s)
  
lib = ffi.dlopen(r"D:\Projects\tomo_cmake\tomo_cmake\build\projector_py\Release\projector_py.dll")

class GeometryFanBeam:

    def __init__(self, angles, i_detectorCount, i_detectorSize, i_dSourceToObj, i_dObjToDetector, imgCenterX, imgCenterY):
        if type(angles) == np.ndarray:
            data = angles.__array_interface__["data"]
            c_angles = ffi.cast("double*", data[0])
        else:
            c_angles = ffi.new("double[]", angles)
        self._handler = ffi.gc(lib.geometry_fan_beam_create(c_angles, len(angles), i_detectorCount, i_detectorSize, i_dSourceToObj,i_dObjToDetector, imgCenterX, imgCenterY), lib.geometry_fan_beam_free)

class GeometryParallel:

    def __init__(self, angles, i_detectorCount, i_detectorSize, imgCenterX, imgCenterY):
        if type(angles) == np.ndarray:
            data = angles.__array_interface__["data"]
            c_angles = ffi.cast("double*", data[0])
        else:
            c_angles = ffi.new("double[]", angles)
        self._handler = ffi.gc(lib.geometry_parallel_create(c_angles, len(angles), i_detectorCount, i_detectorSize, imgCenterX, imgCenterY), lib.geometry_parallel_free)

class Pojector:
    def __init__(self, inputImg, imgSize_x, imgSize_y, geometry, sumAlgorithm):
        if type(inputImg) == np.ndarray:
            inputImg = inputImg.flatten()
            data = inputImg.__array_interface__["data"]
            c_inputImg = ffi.cast("float*", data[0])
        else:
            c_inputImg = ffi.new("float[]", inputImg)
        self._handler = ffi.gc(lib.projector_create(c_inputImg, imgSize_x, imgSize_y, geometry._handler, sumAlgorithm), lib.projector_free)
    
    def get_single_projection(self, angle):
        c_out_arr = ffi.new("float**")
        c_out_size = ffi.new("size_t*")
        lib.projector_get_single_projection(c_out_arr, c_out_size, self._handler, angle)
        data = ffi.gc(c_out_arr[0], lib.float_free)
        return np.array(list(data[0: c_out_size[0]]))
    
    def get_full_projection(self):
        c_out_arr = ffi.new("float**")
        c_out_size_x = ffi.new("size_t*")
        c_out_size_y = ffi.new("size_t*")
        lib.projector_get_full_projection(c_out_arr, c_out_size_x, c_out_size_y, self._handler)
        data = ffi.gc(c_out_arr[0], lib.float_free)
        return np.array(list(data[0: c_out_size_x[0]*c_out_size_y[0]])).reshape((c_out_size_y[0], c_out_size_x[0]))




angles = np.arange(10, dtype = np.float32)/180.0*np.pi
size = 20;
geom = GeometryFanBeam(angles, 50, 1, 7000, 500, size*0.5, size*0.5)
#inputImg = list(range(16))
inputImg = np.arange(size*size, dtype = np.float32).reshape((size,size))*(1.0/(size*size))
proj = Pojector(inputImg, size, size, geom, 0)
a = proj.get_full_projection()
print(a)