__author__ = 'embrun'
import glob
import Image3D
import numpy as np
import SimpleITK as sitk
from math import pi


def command_iteration(method):
   print (method.GetMetricValue())



if __name__ == "__main__":
    print 'K Edge Substraction'

    aboveFolder='/Users/embrun/Experiences/Wiart_MD1011/aboveBin4'
    belowFolder='/Users/embrun/Experiences/Wiart_MD1011/BelowBin4'

    #M1=water
    #M2=Au
    #E1=79.7keV
    #E2=81.7keV

    muM1E1=0.184
    muM1E2=0.182

    muM2E1=2.204
    muM2E2=8.628

    mini=-0.25
    maxi=1.00


    aboveImage=Image3D.Image3D(folderName=aboveFolder)
    aboveImage.createListOfFiles('tif')
    belowImage=Image3D.Image3D(folderName=belowFolder)
    belowImage.createListOfFiles('tif')

    print 'Loading above Image...'
    aboveImage.loadSlices()
    print 'dataType'
    print aboveImage.data.dtype
    print 'Loading below Image...'
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



    print 'Resampling'
    aboveImage.resampleIn32Bit(mini,maxi)
    belowImage.resampleIn32Bit(mini,maxi)

    newNbSlices,newWidth,newHeight=belowImage.getDimensions()
    substractedImage=Image3D.Image3D(nbSlices=newNbSlices,width=newWidth,height=newHeight)

    substractedImage.data= np.subtract(aboveImage.data,belowImage.data)
    #substractedImage.data[substractedImage.data<0.]=0.
    substractedImage.save3DImage('/Users/embrun/Experiences/Wiart_MD1011/WithoutRegistration/brut')


    print 'imResize'
    #aboveImageBin4=Image3D.Image3D(nbSlices=newNbSlices,width=newWidth,height=newHeight)
    #aboveImageBin4.data=np.copy(aboveImage.data)
    #belowImageBin4=Image3D.Image3D(nbSlices=newNbSlices,width=newWidth,height=newHeight)
    #belowImageBin4.data=np.copy(belowImage.data)

    #aboveImageBin4.imresize(0.25)
    #belowImageBin4.imresize(0.25)

    belowImage.save3DImage('/Users/embrun/Experiences/Wiart_MD1011/resampled/below')

    #aboveData=np.copy(aboveImageBin4.data)
    #belowData=np.copy(belowImageBin4.data)


    belowImage.registerImages(aboveImage)


    print 'Conversion to ITK data type'
    #aboveDataITK=sitk.GetImageFromArray(aboveData)
    #belowDataITK=sitk.GetImageFromArray(belowData)

    #belowImageITK=sitk.GetImageFromArray(belowImage.data)

    #print init_t
    #reg_method= sitk.ImageRegistrationMethod()
    #reg_method.SetMetricAsMeanSquares();
    #reg_method.SetOptimizerAsRegularStepGradientDescent( 0.4,0.00001,200,0.5)


    #tx=sitk.TranslationTransform( aboveDataITK.GetDimension() )
    #reg_method.SetInitialTransform(tx )
    #reg_method.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(reg_method) )

    #interpolator=sitk.sitkLinear
    #reg_method.SetInterpolator( interpolator )
    #outTx = reg_method.Execute(aboveDataITK, belowDataITK)
    #print outTx
    #belowImageITKRegistered= sitk.Resample(belowImageITK, aboveDataITK, outTx, interpolator, 0.0,belowImageITK.GetPixelIDValue())
    #belowImageITK=belowImageITKRegistered


    #------------------------------------------
    # reg_method= sitk.ImageRegistrationMethod()
    # reg_method.SetMetricAsMeanSquares();
    # #reg_method.SetOptimizerAsRegularStepGradientDescent( 0.4,0.00001,200,0.5)
    # sample_per_axis=1200
    # reg_method.SetOptimizerAsExhaustive([sample_per_axis//2,0,0])
    # reg_method.SetOptimizerScales([2.0*pi/sample_per_axis, 1.0,1.0])
    #
    # tx=sitk.VersorRigid3DTransform(  )
    # reg_method.SetInitialTransform(tx )
    # reg_method.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(reg_method) )
    #
    # interpolator=sitk.sitkLinear
    # reg_method.SetInterpolator( interpolator )
    # outTx = reg_method.Execute(aboveDataITK, belowDataITK)
    # print outTx
    # belowImageITKRegistered= sitk.Resample(belowImageITK, belowImageITK, outTx, interpolator, 0.0,belowImageITK.GetPixelIDValue())


    #
    # reg_method= sitk.ImageRegistrationMethod()
    # init_t=sitk.CenteredTransformInitializer(aboveDataITK, belowDataITK, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    # print init_t
    # reg_method.SetInitialTransform(init_t)
    # reg_method.SetOptimizerAsRegularStepGradientDescent(0.3,0.001,20,0.6)
    # reg_method.SetMetricAsMeanSquares()
    # reg_method.SetOptimizerScales([1,0.001,0.001,0.001,0.001])
    # outTx = reg_method.Execute(aboveDataITK, belowDataITK)
    # reg_method.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(reg_method) )
    # print outTx
    # belowImageITKRegistered= sitk.Resample(belowImageITK, belowImageITK, outTx, interpolator, 0.0,belowImageITK.GetPixelIDValue())
    #------------------------------------------


    print 'Resampling'
    print 'Converting to numpy array'
    #belowImage.data = sitk.GetArrayFromImage(belowImageITKRegistered)


    print 'Registration'
    #aboveImage.registerWithRef(belowImage)
    print 'Registration Done'
    #aboveImage.save3DImage('/Users/embrun/Essai/RegisteredImage/above')
    #belowImage.save3DImage('/Users/embrun/Essai/RegisteredImage/below')

    #substractedImage=Image3D.Image3D(nbSlices=maxNbSlices,width=maxWidth,height=maxHeight)
    substractedImage.data= np.subtract(aboveImage.data,belowImage.data)
    #substractedImage.data[substractedImage.data<0.]=0.

    # denominateur=muM1E1*muM2E2-muM2E1*muM1E2
    # leftPart=-muM2E1*belowImage.data
    # rightPart=muM2E2*aboveImage.data
    # substractedImage.data= np.add(leftPart,rightPart)
    # substractedImage.divide(denominateur)

    substractedImage.save3DImage('/Users/embrun/Experiences/Wiart_MD1011/TestRegistration/test')









