# coding: utf-8
import numpy as np
from PIL import Image
import math
from tomopy.misc.phantom import shepp2d

def create_phantom(dim, model, size_x, size_y, size_z):
    phantom_path = "phantoms/"+model+"_"+dim
    phantom_bmp_path = phantom_path + ".bmp"
    phantom_txt_path = phantom_path + ".txt"
    phantom = np.zeros((size_x, size_y, size_z))
    
    if model == "phantom":
        phantom2D = np.squeeze(shepp2d(size_x))
    else:
        phantom2D = np.zeros((size_x, size_y))
        phantom2D[5:6, 5:6] = 200
        phantom2D[-6:-5, 5:6] = 150
        phantom2D[5:6, -6:-5] = 100
        phantom2D[-6:-5, -6:-5] = 50
        phantom2D[size_x//2:size_x//2+2, -6:-5] = 50
        phantom2D[-6:-5, size_y//2:size_y//2+2] = 50
        phantom2D[size_x//2:size_x//2+2, size_y//2:size_y//2+2] = 50
        phantom2D[5:6, size_y//2:size_y//2+2] = 50
        phantom2D[size_x//2:size_x//2+2, 5:6] = 50
    
    for i in range(size_z):
        phantom[:, :, i] = phantom2D
    
    with open(phantom_txt_path, "w") as f:
        f.write(str(size_x) + " " + str(size_y) + " " + str(size_z)+"\n")
        for k in range(size_z):
            for j in range(size_y):
                for i in range(size_x):
                    f.write(str(phantom[i, j, k]) + " ")
        f.write("\n")
    print(phantom.shape)
    if model == "phantom":
        phantom = np.floor(phantom2D/phantom2D.max()*255.0)
    else:
        phantom = np.floor(phantom2D)
    im = Image.fromarray(phantom.astype(np.uint8))
    im.convert("L")
    im.save(phantom_bmp_path)
    
def main():
    size_x, size_y, size_z = (512, 512, 100);
    create_phantom("2D", "phantom", size_x, size_y, 1)
    size_x, size_y, size_z = (128, 128, 128);
    create_phantom("3D", "phantom", size_x, size_y, size_z)
    size_x, size_y, size_z = (200, 100, 11);
    create_phantom("2D", "model", size_x, size_y, 1)
    create_phantom("3D", "model", size_x, size_y, size_z)
    
main()
