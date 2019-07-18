'''

This class reads images and saves Transforms

Created 09/01/2018 by Luca Fardin.
Based on the code by L. Broche

'''

import glob
import numpy as np
import fabio
import os
from PIL import Image
from scipy import misc


class InputOutput():

    def __init__(self):
        print('Init I/O')

    def ReadEdfVolume(self, path):

        print("Loading Images ...")
        image_names = glob.glob(path + "/*.edf")
        image_names = sorted(image_names)
        imageobj = fabio.open(image_names[0])
        image = imageobj.data
        image_vol = np.zeros((image.shape[0], image.shape[1], len(image_names)))

        z = 0
        for name_image in image_names:
            print(name_image)
            imageobj = fabio.open(name_image)
            image_vol[:, :, z] = imageobj.data
            z += 1
        # if z==2: break
        return image_vol

    def ReadTifVolume(self, path):

        print("Loading Images ...")
        image_names = glob.glob(path + "/*.tif")
        image_names = sorted(image_names)
        imageobj = Image.open(image_names[0])
        image = np.array(imageobj)
        image_vol = np.zeros((image.shape[0], image.shape[1], len(image_names)))

        z = 0
        for name_image in image_names:
            print(name_image)
            imageobj = Image.open(name_image)
            image_vol[:, :, z] = np.array(imageobj)
            z += 1
        # if z==2: break
        return image_vol

    def SaveEdf(self, output_image, path, root):
        print("Saving Images ...")

        # Check Folder Exists
        if not os.path.exists(path):
            os.makedirs(path)

        for i in range(output_image.shape[2]):
            filename = path + root + str(i).zfill(4) + '.edf'
            fabio.edfimage.EdfImage(data=output_image[:, :, i].astype(np.float32), header={}).save(filename)

    def ReadPngVolume(self, path):

        print("Loading Images ...")
        image_names = glob.glob(path + "/*.png")
        image_names = sorted(image_names)
        image = misc.imread(image_names[0])
        print(image.shape)
        image_vol = np.zeros((image.shape[0], image.shape[1], len(image_names)))

        z = 0
        for name_image in image_names:
            print(name_image)
            image = misc.imread(name_image)
            image_vol[:, :, z] = image
            z += 1
        # if z==2: break
        return image_vol.astype(np.uint8)

    def SavePng(self, output_image, path, root):
        print("Saving Images ...")

        # Check Folder Exists
        if not os.path.exists(path):
            os.makedirs(path)

        for i in range(output_image.shape[2]):
            filename = path + root + str(i).zfill(4) + '.png'
            misc.imsave(filename, output_image[:, :, i])

    def SaveTiff(self, output_image, path, root):
        print("Saving Images ...")

        # Check Folder Exists
        if not os.path.exists(path):
            os.makedirs(path)

        for i in range(output_image.shape[2]):
            filename = path + root + str(i).zfill(4) + '.tif'
            misc.imsave(filename, output_image[:, :, i])

    def ReadEdfTest(self, path):

        print("Loading Images ...")
        image_names = glob.glob(path + "/*.edf")
        image_names = sorted(image_names)
        imageobj = fabio.open(image_names[0])
        imageobj = imageobj.data
        image_vol = np.zeros((imageobj.shape[0], imageobj.shape[1], 10))

        z = 0
        for name_image in image_names:
            print(name_image)
            imageobj = fabio.open(name_image)
            image_vol[:, :, z] = imageobj.data
            z += 1
            if z == 10: break
        return image_vol

    def ReadTifTest(self, path):

        print("Loading Images ...")
        image_names = glob.glob(path + "/*.tif")
        image_names = sorted(image_names)
        imageobj = Image.open(image_names[0])
        image = np.array(imageobj)
        image_vol = np.zeros((image.shape[0], image.shape[1], 5))

        z = 0
        for name_image in image_names:
            print(name_image)
            imageobj = Image.open(name_image)
            image_vol[:, :, z] = np.array(imageobj)
            z += 1
            if z == 5: break
        return image_vol
