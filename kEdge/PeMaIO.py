

import EdfFile as edf
from PIL import Image
import numpy as np
from scipy.ndimage import fourier_shift
import glob
import array
import scipy


def readImage(filename):
    if filename.endswith('.edf') or filename.endswith('.tiff') :
        edfFile = edf.EdfFile(filename, access='rb')
        im2D = edfFile.GetData(0)

    if filename.endswith('.tiff') or filename.endswith('.tif') or filename.endswith('.png') or  filename.endswith('.TIF')  or filename.endswith('.TIFF') :
        img=Image.open(filename)
        im2D= np.array(img)

    return im2D


def writeImage(filename, data):
    if filename.endswith('.edf') or filename.endswith('.tiff') :
        edfFile = edf.EdfFile(filename, access='wb-')
        edfFile.WriteImage({}, data)

    else:
        typeImage=data.dtype
        if typeImage==np.uint8 :
            scipy.misc.imsave(filename, data)

        if typeImage==np.bool :
            toStore=np.zeros(data.shape,dtype=np.uint8)
            toStore[data==True]=255
            scipy.misc.imsave(filename, toStore)

        if typeImage==np.float32 or typeImage==np.float16 or typeImage==np.float64 :
            scipy.misc.imsave(filename, data)

    return im2D
