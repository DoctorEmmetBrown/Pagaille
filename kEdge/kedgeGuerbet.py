__author__ = 'embrun'
from skimage.feature import register_translation
import EdfFile as edf
import glob
import string
import os
import numpy as np
import SimpleITK as sitk
import Image3D








if __name__ == "__main__":
    print 'On demarre'
    belowFolder='/Volumes/BM05/imaging/embrun/160307_gado/ReconstructionPhantom/below/'
    aboveFolder='/Volumes/BM05/imaging/embrun/160307_gado/ReconstructionPhantom/above/'
    outputFolder='/Volumes/BM05/imaging/embrun/160307_gado/ReconstructionPhantom/Substracted/'


    muM1E1=0.250
    muM1E2=0.247

    muM2E1=30.782
    muM2E2=142.93

    aboveImage=Image3D.Image3D(folderName=aboveFolder)
    aboveImage.createListOfFiles('edf')
    belowImage=Image3D.Image3D(folderName=belowFolder)
    belowImage.createListOfFiles('edf')

    aboveImage.loadSlices()
    belowImage.loadSlices()

    aboveNbSlices,aboveWidth,aboveHeight=aboveImage.getDimensions()
    belowNbSlices,belowWidth,belowHeight=belowImage.getDimensions()

    substractedImage=Image3D.Image3D(nbSlices=aboveNbSlices,width=belowWidth,height=belowHeight)

    denominateur=muM1E1*muM2E2-muM2E1*muM1E2
    leftPart=-muM1E1*belowImage.data
    rightPart=muM1E1*aboveImage.data
    substractedImage.data= np.add(leftPart,rightPart)
    substractedImage.divide(denominateur)

    substractedImage.save3DImage(outputFolder+'/substracted_')










