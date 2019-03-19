__author__ = 'embrun'

import EdfFile as edf
from PeMaIO import readImage as readIm
import numpy as np
from scipy.ndimage import fourier_shift
import glob
import array
#from skimage.feature import register_translation
import scipy
from skimage.transform import resize
import math
import SimpleITK as sitk


def command_iteration(method):
   print (method.GetMetricValue())


class Image3D:
    def __init__(self, *args, **kwargs):
        if 'scanNumber' in kwargs:
            self.scanNumber = kwargs.pop('scanNumber')
        else:
            self.scanNumber = 0

        if 'nbSlices' in kwargs:
            self.nbSlices = kwargs.pop('nbSlices')
        else:
            self.nbSlices = -1

        if 'folderName' in kwargs:
            self.folderName = kwargs.pop('folderName')
        else:
            self.folderName = None

        if 'listofFiles' in kwargs:
            self.listofFiles = kwargs.pop('listofFiles')
        else:
            self.listofFiles = []

        if 'dtype' in kwargs:
            self.dtype = kwargs.pop('dtype')
        else:
            self.dtype = np.float32

        self.data = None

        self.width = -1
        self.height = -1
        if 'width' in kwargs:
            self.width = kwargs.pop('width')
        if 'height' in kwargs:
            self.height = kwargs.pop('height')

            # if self.width>0:
            #    self.createImage()

    def isEmptyScan(self):
        if not self.listofFiles:
            return True
        else:
            return False

    def createListOfFiles(self, extension):
        self.listofFiles = glob.glob(self.folderName + '/*' + extension)
        self.listofFiles.sort()
        self.nbSlices = len(self.listofFiles)

    def appendFileName(self, num, ff):
        self.listofFiles.append(ff)

    def printListOfFiles(self):
        if not self.listofFiles:
            print('Empty List')
        for ff in self.listofFiles:
            print(ff)

    def createImage(self):
        self.data = np.zeros((self.nbSlices, self.width, self.height), self.dtype)

    def setDimensions(self):
        firstFile = self.listofFiles[0]
        print(firstFile)
        im2D = readIm(firstFile)
        self.width, self.height = im2D.shape
        self.dtype = im2D.dtype
        self.data = np.zeros((self.nbSlices, self.width, self.height), im2D.dtype)

    def resampleIn32Bit(self,mini,maxi):
        #newdata= np.zeros((self.nbSlices, self.width, self.height), np.float32)
        maxImage=np.amax(self.data)
        minImage= np.amin(self.data)
        #print('RESAMPLE dans image:'+str(minImage)+' maxImage:'+str(maxImage))
        #print('Conversion : '+str(mini)+' maxi:'+str(maxi))
        self.data=np.asarray(self.data,np.float32)
        self.data= (maxi-mini)*(self.data/maxImage) + mini #on a remplacer 65535 par maxImage
        maxImage = np.amax(self.data)
        minImage = np.amin(self.data)
        #print('Apres Resampling dans image:' + str(minImage) + ' maxImage:' + str(maxImage))

        #=newdata


    def loadSlices(self):
        """

        Returns
        -------
        object
        """
        print('Load ')
        if (self.width < 0):
            self.setDimensions()

        i = 0
        for fileN in self.listofFiles:
            fileName = fileN
            #edfFile = edf.EdfFile(fileName, access='rb')
            im2D = readIm(fileName)
            self.data[i, :, :] = im2D
            i += 1

    def sumImage(self, other):
        self.data = np.add(self.data, other.data)

    def divide(self, constant):
        self.data = np.true_divide(self.data, constant)

    def registerWithRef(self, imRef):

        shift, error, diffphase = register_translation(imRef.data, self.data)
      
        self.data = np.roll(self.data, shift=int(shift[0]), axis=0)
        self.data = np.roll(self.data, shift=int(shift[1]), axis=1)
        self.data = np.roll(self.data, shift=int(shift[2]), axis=2)

    def save3DImage(self, outputName):
        self.nbSlices, self.width, self.height = self.data.shape
        for k in range(self.nbSlices):
            textNumSlice = '%4.4d' % k
            finalOutputName = outputName + textNumSlice + '.edf'
            filetoWrite = edf.EdfFile(finalOutputName, access='wb+')
            dataToStore = self.data[k, :, :].squeeze()
            filetoWrite.WriteImage({}, dataToStore)

    def saveSino(self, outputName):
        

        for i in range(self.width):
            textNumSlice = '%4.4d' % i
            finalOutputName = outputName + textNumSlice + '.edf'
            filetoWrite = edf.EdfFile(finalOutputName, access='wb+')
            dataToStore = self.data[:, i, :].squeeze()
            filetoWrite.WriteImage({}, dataToStore)

    def getDimensions(self):
        return self.nbSlices, self.width, self.height

    def redimension(self, newNbSlices, newWidth, newHeight):
        temp = np.zeros((newNbSlices, newWidth, newHeight), self.dtype)
        temp[0:self.nbSlices, 0:self.width, 0:self.height] = self.data
        self.data = temp
        self.nbSlices = newNbSlices
        self.width = newWidth
        self.height = newHeight

    def stackRegisterWithRef(self, imRef):
       
        for k in range(self.nbSlices):
            slice = self.data[k, :, :].squeeze()
            sliceRef = imRef.data[k, :, :].squeeze()
            shift, error, diffphase = register_translation(slice, sliceRef, 100)
            print(shift)

    def zeros(self, newNbSlices, newWidth, newHeight):
        self.data = np.zeros((newNbSlices, newWidth, newHeight), self.dtype)
        self.nbSlices = newNbSlices
        self.width = newWidth
        self.height = newHeight


    def imresize(self, factor):
        self.nbSlices=int(self.nbSlices*factor)
        self.width=int(self.width*factor)
        self.height=int(self.height*factor)
        dim=(self.nbSlices,self.width,self.height)
        self.data=resize(self.data,dim)

    def registerImages(self,fixedImage):
        print('Conversion in ITK format')
        fixedImageITK=sitk.GetImageFromArray(fixedImage.data)
        movingImageITK=sitk.GetImageFromArray(self.data)
        print('Conversion made')
        reg_method= sitk.ImageRegistrationMethod()
        reg_method.SetMetricAsMeanSquares();
        reg_method.SetOptimizerAsRegularStepGradientDescent( 0.4,0.00001,200,0.5)
        tx=sitk.TranslationTransform(fixedImageITK.GetDimension() )
        reg_method.SetInitialTransform(tx )
        reg_method.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(reg_method) )

        interpolator=sitk.sitkLinear
        reg_method.SetInterpolator( interpolator )
        outTx = reg_method.Execute(fixedImageITK, movingImageITK)
        print(outTx)
        movedImageITK= sitk.Resample(movingImageITK, fixedImageITK, outTx, interpolator, 0.0,fixedImageITK.GetPixelIDValue())
        self.data = sitk.GetArrayFromImage(movedImageITK)



    def rigidRegisterImages(self,fixedImage):
        print('Conversion in ITK format')
        fixedImageITK=sitk.GetImageFromArray(fixedImage.data)
        movingImageITK=sitk.GetImageFromArray(self.data)
        print('Conversion made')
        reg_method= sitk.ImageRegistrationMethod()
        init_t=sitk.CenteredTransformInitializer(fixedImageITK, movingImageITK, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
        print('Init transformation')
        print(init_t)
        reg_method.SetInitialTransform(init_t)
        reg_method.SetMetricAsMeanSquares();
        reg_method.SetOptimizerAsRegularStepGradientDescent( 0.4,0.00001,200,0.5)
        reg_method.SetOptimizerScales([1,0.001,0.001,0.001,0.001])
        reg_method.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(reg_method))
        interpolator=sitk.sitkLinear
        reg_method.SetInterpolator( interpolator )
        outTx = reg_method.Execute(fixedImageITK, movingImageITK)
        print(outTx)
        movedImageITK= sitk.Resample(movingImageITK, fixedImageITK, outTx, interpolator, 0.0,fixedImageITK.GetPixelIDValue())
        self.data = sitk.GetArrayFromImage(movedImageITK)




