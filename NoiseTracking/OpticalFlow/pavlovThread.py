from threading import Thread
from NoiseTracking.OpticalFlow import PavlovSingleMaterial
import os
from InputOutput.pagailleIO import openImage,saveEdf
import numpy as np




class OpticalFlowSolverThread(Thread):
    def __init__(self,listOfProjections,listOfSampleFileNames,referenceName,outputFolder):
        Thread.__init__(self)
        self.referenceName = referenceName
        self.listOfProjection=listOfProjections
        self.output=outputFolder
        self.listOfSampleFileNames=listOfSampleFileNames


    def run(self):
        beta = 2.72809e-12
        # beta=1.48E-15
        # n=1-delta +i*beta
        mu = 1.4805e0  # *30.

        # mu=2*k*beta
        # by a constant to produce a low frequency filter
        # thickness_LFF_x30 - for mu *30
        # thickness_LFF_x20 - for mu *20
        # thickness_LFF_x10 - for mu *10
        # thickness_LFF_x5 - for mu *5
        delta = 1.04842e-7
        # delta = 1.32e-7
        # /39.4784176;
        z1 = 0.9  # (m) distance for the msk to the object
        z2 = 2.0  # (m) distance from the object to the detector
        R1 = 42.  # (m) ~ distance from the source to the satellite building from ID17 website
        R2 = 2.0  # (m) Propagation distance from the exit surface of the object
        pix_size = 6.1e-6  # 6 um as in Paganin et al paper

        E = 52
        low_frequency_filter = 1.
        scale = 20.

        Ir = openImage(self.referenceName)
        Ir = np.asarray(Ir, np.float32)

        projectionFiles = self.listOfSampleFileNames
        for numeroProjection in self.listOfProjection:
            numeroProjection = int(numeroProjection)
            print('Processing ' + str(numeroProjection))
            projectionFileName=projectionFiles[numeroProjection]
            Is=openImage(projectionFileName)

            Image_old = np.true_divide(Is, Ir)
            Image_old = 1 - Image_old
            # New smaller array containing image of fibres only

            Image_new = Image_old

            # Calculation of the average value of the image
            average_image = np.mean(Image_new)
            # Correction on the average image. Now the average of the new array is ~0
            Image_new = Image_new - average_image
            img_in = Image_new

            bg_val = 0

            img_out = PavlovSingleMaterial.tie_hom_KMP2Last(img_in, E, R1, R2, pix_size, delta, beta * low_frequency_filter, bg_val, scale)
            img_out = img_out * 1e6
            saveEdf(img_out, self.output+ '/' + os.path.basename(projectionFileName))

