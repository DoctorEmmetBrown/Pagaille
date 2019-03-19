__author__ = 'embrun'
import glob
import Image3D
import numpy as np


if __name__ == "__main__":
    print 'K Edge Substraction'
    belowFolder='/Users/embrun/Essai/Below/'
    aboveFolder='/Users/embrun/Essai/Above/'

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

    maxNbSlices=max(aboveNbSlices,belowNbSlices)
    maxWidth=max(aboveWidth,belowWidth)
    maxHeight=max(aboveHeight,belowHeight)

    print 'maxNbSlices:'+str(maxNbSlices)
    print 'maxWidth:'+str(maxWidth)
    print 'maxHeight:'+str(maxHeight)

    print 'above'
    aboveImage.redimension(maxNbSlices,maxWidth,maxHeight)
    print 'below'
    belowImage.redimension(maxNbSlices,maxWidth,maxHeight)


    print 'Registration'
    aboveImage.registerWithRef(belowImage)
    print 'Registration Done'
    aboveImage.save3DImage('/Users/embrun/Essai/RegisteredImage/above')
    belowImage.save3DImage('/Users/embrun/Essai/RegisteredImage/below')

    substractedImage=Image3D.Image3D(nbSlices=maxNbSlices,width=maxWidth,height=maxHeight)

    denominateur=muM1E1*muM2E2-muM2E1*muM1E2
    leftPart=-muM1E1*belowImage.data
    rightPart=muM1E1*aboveImage.data
    substractedImage.data= np.add(leftPart,rightPart)
    substractedImage.divide(denominateur)

    substractedImage.save3DImage('/Users/embrun/Essai/kedgeSubstracted/ConcentrationGado')









