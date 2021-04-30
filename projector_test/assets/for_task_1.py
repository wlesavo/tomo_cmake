# coding: utf-8
import numpy as np
import pylab as plt
from PIL import Image
import astra
import math
#import tomopy
from tomopy.misc.phantom import shepp2d

def fp(pg, vg, v, algo):

    proj_id = astra.create_projector(algo, pg, vg);
    v_id = astra.data2d.create('-vol', vg, v)
    rt_id = astra.data2d.create('-sino', pg)
    fp_cfg = astra.astra_dict('FP')
    fp_cfg['ProjectorId'] = proj_id
    fp_cfg['VolumeDataId'] = v_id
    fp_cfg['ProjectionDataId'] = rt_id
    fp_id = astra.algorithm.create(fp_cfg)
    astra.algorithm.run(fp_id)
    out = astra.data2d.get(rt_id)
    astra.algorithm.delete(fp_id)
    astra.data2d.delete(rt_id)
    astra.data2d.delete(v_id)
    return out

def main(geometry, algo, model):
    
    phantom_path = "phantom_"+model+"_"+geometry+"_"+algo+".bmp"
    sinogram_path = "sinogram_"+model+"_"+geometry+"_"+algo+".bmp"
    sinogram_txt_path = "sinogram_"+model+"_"+geometry+"_"+algo+".txt"

    size = 256
    angles = np.linspace(0, 2*np.pi, 360, False)
    print(angles[3])
    source_object = 7000
    object_det = 500
    phantom = np.zeros((size, size))
    phantom[5:6, 5:6] = 200
    phantom[-6:-5, 5:6] = 200
    phantom[5:6, -6:-5] = 200
    phantom[-6:-5, -6:-5] = 200
    #for i in range(3, size - 3):
    #    phantom[i,i] = 100
    #    phantom[size - 1 - i,i] = 100
    phantom[size//2:size//2+2, -6:-5] = 50
    phantom[-6:-5, size//2:size//2+2] = 50
    phantom[size//2:size//2+2, size//2:size//2+2] = 50
    phantom[5:6, size//2:size//2+2] = 50
    phantom[size//2:size//2+2, 5:6] = 50
    vol_geom = astra.create_vol_geom((size, size))
    if model=="phantom":
        phantom = np.squeeze(shepp2d(size))
    
    if geometry == "par":
        proj_geom = astra.create_proj_geom("parallel", 1.0, size, angles)
    elif geometry == "fan":
        proj_geom = astra.create_proj_geom("fanflat", 1.0, size, angles, source_object, object_det)
    
    if algo == "area" and geometry == "par":
        algorithm = 'strip'
    elif algo == "area" and geometry == "fan":
        algorithm = "strip_fanflat"
    elif algo == "line" and geometry == "par":
        algorithm = 'line'
    elif algo == "line" and geometry == "fan":
        algorithm = 'line_fanflat'
        
    sigogram = fp(proj_geom, vol_geom, phantom, algorithm)
    
    with open("angles.txt", "w") as f :
        f.write(str(angles.size))
        f.write("\n")
        for i in range(len(angles)):
            f.write(str(angles[i]))
            f.write("\n")
    with open(sinogram_txt_path, "w") as f :
        for i in range(len(angles)):
            f.write(" ".join([str(x) for x in sigogram[i, :].tolist()]))
            f.write("\n")
    print(sigogram.max())
    sigogram = np.floor(sigogram/sigogram.max()*255.0)
    #sigogram = np.transpose(sigogram)
    #print(sigogram[0,0], sigogram[359,0], sigogram[0,255], sigogram[359, 255])
    im = Image.fromarray(phantom.astype(np.uint8))
    im.convert("L")
    im.save(phantom_path)
    im = Image.fromarray(np.transpose(sigogram.astype(np.uint8)))
    im.convert("L")
    im.save(sinogram_path)

def create_all():
    main("par", "line", "phantom")
    main("fan", "line", "phantom")
    main("par", "area", "phantom")
    main("fan", "area", "phantom")
    main("par", "line", "model")
    main("fan", "line", "model")
    main("par", "area", "model")
    main("fan", "area", "model")
    
create_all()
