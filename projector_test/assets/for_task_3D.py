# coding: utf-8
import numpy as np
from PIL import Image
import astra
import math
#import tomopy

def fp(proj_geom, vol_geom, volume):

    sinogram_id = astra.data3d.create('-sino', proj_geom)
    # Create a data object for the reconstruction
    rec_id = astra.data3d.create('-vol', vol_geom, data=volume)

    # Set up the parameters for a reconstruction algorithm using the GPU
    cfg = astra.astra_dict('FP3D_CUDA')
    cfg['VolumeDataId'] = rec_id
    cfg['ProjectionDataId'] = sinogram_id
    cfg['option'] = {}

    alg_id = astra.algorithm.create(cfg)
    astra.algorithm.run(alg_id, 1)
    res_sino = astra.data3d.get(sinogram_id)
    
    astra.algorithm.delete(alg_id)
    astra.data3d.delete(rec_id)
    astra.data3d.delete(sinogram_id)
    astra.clear()
    return res_sino

def main():
    
    phantom_path = "shepp2D"+".bmp"
    sinogram_path = "sinogram_3D"+".bmp"
    phantom_txt_path = "phantom_3D"+".txt"
    sinogram_txt_path = "sinogram_3D"+".txt"
    
    size_x, size_y, size_z = (256, 256, 20);
    detector_size_x, detector_size_y = (256, 40)
    #detector_size_x, detector_size_y = (40, 30)
    angles = np.linspace(0, np.pi, 10, False)
    
    source_object = 200
    object_det = 200
    im = Image.open(phantom_path)
    phantom2D = np.zeros((size_x, size_y))
    pixels = im.load()
    for j in range(size_y):
        for i in range(size_x):
            phantom2D[i, j] = pixels[i, j]
    #size_x, size_y, size_z = (100, 50, 20);
    phantom = np.zeros((size_z, size_y, size_x, ))
    for k in range(size_z):
        phantom[k,:,:] = phantom2D*abs(size_z - (size_z//2 - k))    
    #phantom[5:6, 5:6, 5:6]       = 10
    #phantom[-6:-5, 5:6, 5:6]   = 20
    #phantom[5:6, -6:-5, 5:6]     = 30
    #phantom[-6:-5, -6:-5, 5:6] = 40
    #phantom[5:6, 5:6, -6:-5]       = 50
    #phantom[-6:-5, 5:6, -6:-5]   = 60
    #phantom[5:6, -6:-5, -6:-5]     = 70
    #phantom[-6:-5, -6:-5, -6:-5] = 80
    #for k in range(size_z):
    #    for j in range(size_x):
    #        for i in range(size_y):
    #            if i>5 and i < 10:
    #                phantom[k, i, j] = 100
    #            if j>5 and j < 10:
    #                phantom[k, i, j] = 100
                
                #if (j+k+i)%5 < 2:
                #   phantom[k, i, j] += 10
    with open(phantom_txt_path, "w") as f :
        f.write(str(size_x) + " " + str(size_y) + " " + str(size_z)+"\n")
        for k in range(size_z):
            for j in range(size_y):
                for i in range(size_x):
                    f.write(str(phantom[size_z - k - 1, j, size_x - i - 1]) + " ")
        f.write("\n")

    vol_geom = astra.create_vol_geom((size_y, size_x, size_z))
    proj_geom = astra.create_proj_geom('cone',  1.0, 1.0, detector_size_y, detector_size_x, angles, source_object, object_det);
    sinogram = fp(proj_geom, vol_geom, phantom)
    sino_arr = np.zeros((detector_size_x, angles.shape[0]*detector_size_y))
    det_count = detector_size_y*detector_size_x
    with open(sinogram_txt_path, "w") as f :
        f.write(str(angles.shape[0]) + " " + str(detector_size_x) + " " + str(detector_size_y)+"\n")
        for k in range(angles.shape[0]):
            for j in range(detector_size_y):
                for i in range(detector_size_x):
                    sino_arr[i, j + k * detector_size_y] = sinogram[j, k, i]
                    f.write(str(sinogram[j, k, i]) + " ")
        f.write("\n")
    
    
    
    print(sinogram.max())
    
    print(sinogram.shape)
    #sinogram = np.transpose(sinogram)
    #print(sinogram[0,0], sinogram[359,0], sinogram[0,255], sinogram[359, 255])
    sino_arr = np.floor(sino_arr.T/sino_arr.max()*255.0)
    im = Image.fromarray(sino_arr.astype(np.uint8))
    im.convert("L")
    im.save(sinogram_path)
  
main()
