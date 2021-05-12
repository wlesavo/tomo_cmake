# coding: utf-8
import numpy as np
from PIL import Image
import astra
import math
from timeit import default_timer as tmr

def fp(pg, vg, v, algo):
    proj_id = astra.create_projector(algo, pg, vg);
    v_id = astra.data2d.create('-vol', vg, v)
    rt_id = astra.data2d.create('-sino', pg)
    fp_cfg = astra.astra_dict('FP')
    fp_cfg['ProjectorId'] = proj_id
    fp_cfg['VolumeDataId'] = v_id
    fp_cfg['ProjectionDataId'] = rt_id
    fp_id = astra.algorithm.create(fp_cfg)
    start = tmr()
    astra.algorithm.run(fp_id)
    end = tmr()
    print("time ms: ", (end - start)*1000, "\n")
    out = astra.data2d.get(rt_id)
    astra.algorithm.delete(fp_id)
    astra.data2d.delete(rt_id)
    astra.data2d.delete(v_id)
    return out

def fp_3D(proj_geom, vol_geom, volume):

    sinogram_id = astra.data3d.create('-sino', proj_geom)
    # Create a data object for the reconstruction
    rec_id = astra.data3d.create('-vol', vol_geom, data=volume)

    # Set up the parameters for a reconstruction algorithm using the GPU
    cfg = astra.astra_dict('FP3D_CUDA')
    cfg['VolumeDataId'] = rec_id
    cfg['ProjectionDataId'] = sinogram_id
    cfg['option'] = {}

    alg_id = astra.algorithm.create(cfg)
    start = tmr()
    astra.algorithm.run(alg_id, 1)
    end = tmr()
    print("time ms: ", (end - start)*1000, "\n")
    res_sino = astra.data3d.get(sinogram_id)
    
    astra.algorithm.delete(alg_id)
    astra.data3d.delete(rec_id)
    astra.data3d.delete(sinogram_id)
    astra.clear()
    return res_sino

def calculate_sino_2D(geometry, algo, model):
    
    phantom_path = "phantoms/"+model+"_2D" + ".txt"
    sinogram_bmp_path = "target_sinograms/sinogram_2D_"+model+"_"+geometry+"_"+algo+".bmp"
    sinogram_txt_path = "target_sinograms/sinogram_2D_"+model+"_"+geometry+"_"+algo+".txt"
    
    detector_size_x = 256
    angles = np.linspace(0, 2*np.pi, 360, False)
    source_object = 7000
    object_det = 1000
        
    with open(phantom_path, "r") as f:
        size_x, size_y, size_z = [int(x) for x in f.readline().split()]
        phantom = np.zeros((size_x, size_y))
        a = [float(x) for x in f.readline().split()]
        for j in range(size_y):
            for i in range(size_x):
                coor = i + size_x * j
                phantom[i, j] = a[coor]
    
    print("2D", algo, geometry)
    print("img size:", size_x, size_y)
    print("detector size:", detector_size_x)
    print("angle count:", angles.shape[0])
    
    vol_geom = astra.create_vol_geom((size_x, size_y))
    
    if geometry == "par":
        proj_geom = astra.create_proj_geom("parallel", 1.0, detector_size_x, angles)
    elif geometry == "fan":
        proj_geom = astra.create_proj_geom("fanflat", 1.0, detector_size_x, angles, source_object, object_det)
    
    if algo == "area" and geometry == "par":
        algorithm = 'strip'
    elif algo == "area" and geometry == "fan":
        algorithm = "strip_fanflat"
    elif algo == "line" and geometry == "par":
        algorithm = 'line'
    elif algo == "line" and geometry == "fan":
        algorithm = 'line_fanflat'
        
    sinogram = fp(proj_geom, vol_geom, phantom, algorithm)
    
    with open(sinogram_txt_path, "w") as f :
        for i in range(len(angles)):
            f.write(" ".join([str(x) for x in sinogram[i, :].tolist()]))
            f.write("\n")
    sinogram = np.floor(sinogram/sinogram.max()*255.0)
    im = Image.fromarray(np.transpose(sinogram.astype(np.uint8)))
    im.convert("L")
    im.save(sinogram_bmp_path)

def calculate_sino_3D(model, geometry):
    
    phantom_path = "phantoms/"+model+"_3D"+".txt"
    sinogram_txt_path = "target_sinograms/sinogram_3D_"+model+".txt"
    sinogram_bmp_path = "target_sinograms/sinogram_3D_"+model+".bmp"
    detector_size_x, detector_size_y = (128, 128)
    angles = np.linspace(0, np.pi, 180, False)
    source_object = 7000
    object_det = 1000
    
    with open(phantom_path, "r") as f:
        size_x, size_y, size_z = [int(x) for x in f.readline().split()]
        phantom = np.zeros((size_z, size_y, size_x))
        a = [float(x) for x in f.readline().split()]
        for k in range(size_z):
            for j in range(size_y):
                for i in range(size_x):
                    coor = i + size_x * j + size_x * size_y * k
                    phantom[k, j, i] = a[coor]
    
    print("3D", "line", geometry)
    print("img size:", size_x, size_y, size_z)
    print("angle count:", angles.shape[0])
    print("detector size:", detector_size_x, detector_size_y)
    
    vol_geom = astra.create_vol_geom((size_y, size_x, size_z))
    if geometry == "par":
        proj_geom = astra.create_proj_geom('parallel3d', 1.0, 1.0, detector_size_y, detector_size_x, angles);
    elif geometry == "fan":
        proj_geom = astra.create_proj_geom('cone', 1.0, 1.0, detector_size_y, detector_size_x, angles, source_object, object_det);
    
 
    sinogram = fp_3D(proj_geom, vol_geom, phantom)
    sino_arr = np.zeros((detector_size_x, angles.shape[0]*detector_size_y))
    det_count = detector_size_y*detector_size_x
    with open(sinogram_txt_path, "w") as f :
        f.write(str(angles.shape[0]) + " " + str(detector_size_x) + " " + str(detector_size_y)+"\n")
        for k in range(angles.shape[0]):
            for j in range(detector_size_y):
                for i in range(detector_size_x):
                    sino_arr[i, j + k * detector_size_y] = sinogram[j, k, min(detector_size_x - i - 1, detector_size_x - 1)]
                    f.write(str(sinogram[j, k, min(detector_size_x - i - 1, detector_size_x - 1)]) + " ")
        #f.write("\n")
    sino_arr = np.floor(sino_arr.T/sino_arr.max()*255.0)
    im = Image.fromarray(sino_arr.astype(np.uint8))
    im.convert("L")
    im.save(sinogram_bmp_path)

def create_all():
    calculate_sino_2D("par", "line", "phantom")
    calculate_sino_2D("fan", "line", "phantom")
    calculate_sino_2D("par", "area", "phantom")
    calculate_sino_2D("fan", "area", "phantom")
    calculate_sino_2D("par", "line", "model")
    calculate_sino_2D("fan", "line", "model")
    calculate_sino_2D("par", "area", "model")
    calculate_sino_2D("fan", "area", "model")
    #calculate_sino_3D("model", "fan")
    #calculate_sino_3D("model", "par")
    calculate_sino_3D("phantom", "fan")
    calculate_sino_3D("phantom", "par")
    
    
   
create_all()
#calculate_sino_3D("model")