import Image3D
import numpy as np
from ssim import SSIM

def command_iteration(method):
   print (method.GetMetricValue())



if __name__ == "__main__":
    print 'K Edge Substraction'

    aboveFolder='/Users/embrun/Experiences/Wiart_MD1011/R0842_03_aboveAu/'
    belowFolder='/Users/embrun/Experiences/Wiart_MD1011/R0842_03_belowAu/'

    muM1E1=0.184
    muM1E2=0.182

    muM2E1=2.204
    muM2E2=8.628

    mini=-0.25
    maxi=1.00


    aboveImage=Image3D.Image3D(folderName=aboveFolder)
    aboveImage.createListOfFiles('edf')
    belowImage=Image3D.Image3D(folderName=belowFolder)
    belowImage.createListOfFiles('edf')

    print 'Loading above Image...'
    aboveImage.loadSlices()
    print 'dataType'
    print aboveImage.data.dtype
    print 'Loading below Image...'
    belowImage.loadSlices()

    newNbSlices,newWidth,newHeight=belowImage.getDimensions()
    substractedImage=Image3D.Image3D(nbSlices=newNbSlices,width=newWidth,height=newHeight)

    substractedImage.data= np.subtract(aboveImage.data,belowImage.data)
    substractedImage.save3DImage('/Users/embrun/Experiences/Wiart_MD1011/WithoutRegistration/brut')

    print 'Registration'
    belowImage.registerImages(aboveImage)

    print 'Registration Done'
    substractedImage.data= np.subtract(aboveImage.data,belowImage.data)

    substractedImage.save3DImage('/Users/embrun/Experiences/Wiart_MD1011/TestRegistration/test')









