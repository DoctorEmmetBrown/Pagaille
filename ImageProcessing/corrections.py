
from scipy.ndimage import gaussian_filter
from scipy.ndimage import fourier_shift
import numpy as np
import spytIO
import EdfFile as edf




def registerRefAndSample(Ir,Is,precision):
    from skimage.feature import register_translation  as register_translation
    # correction of beam movementn
    shifts, error, phasediff = register_translation(Is, Ir, precision)
    offset_image = fourier_shift(np.fft.fftn(Ir), shifts)
    offset_image = np.fft.ifftn(offset_image)
    print (' the register correction shift is ( in pixels ) : ')
    print (shifts)
    print (' ----------------------- registration correction done------------------------ ')
    return offset_image.real


def getShift(Im1,Im2,precision=1000):
    shifts, error, phasediff = register_translation(Im1, Im2, precision)
    return shifts




def registerImagesBetweenThemselves(im1, im2,sliceReference=0):
    from skimage.feature import register_translation  as register_translation
    Is=im1.copy()
    Ir= im2.copy()

    refImage=Is[sliceReference,:,:]
    nbSlices,height,width=Is.shape

    for slice in (range(0,nbSlices)):
        if slice != sliceReference:
            imageToMove=Is[slice,:,:]
            shift=getShift(refImage,imageToMove,precision=1000)

            if shift[0]>1:
                shift[0]=0

            print('Slice:'+str(slice)+'Shift:'+str(shift))
            offset_image = fourier_shift(np.fft.fftn(imageToMove), shift)
            offset_image = np.fft.ifftn(offset_image).real
            Is[slice,:,:]=offset_image

            refToMove=Ir[slice,:,:]
            offset_image = fourier_shift(np.fft.fftn(refToMove), shift)
            offset_image = np.fft.ifftn(offset_image).real
            Ir[slice, :, :] = offset_image

    return Is,Ir



def errorDetection (Ir, Is):
    # correction of failed acquired  couple Ir and Is
    # Ir and Is are 3D array as [nbpoint,height,width] because openSeq return toReturn = np.zeros((len(filenames), height, width))
    nbPoints,height,width=Ir.shape
    print(len(nbPoints))
    pointsToDelete=[]
    for coupleOfPopints in range (0,nbPoints):
        sample=Is[coupleOfPopints,:,:]
        reference = Ir[coupleOfPopints, :, :]
        shifts, error, phasediff = registerRefAndSample(reference,sample, 1000)
        if shifts[0]> 10 or shifts [1] >10:
            pointsToDelete.append(coupleOfPopints)
            print('ERROR couple : deletion in progress ')

    if (len(pointsToDelete)>0):
        Ir=np.delete(Ir, pointsToDelete,axis=0)
        Is = np.delete(Is, pointsToDelete, axis=0)

    return Ir,Is



    print(' ----------------------- error detection correction done------------------------ ')


def normalization (Im, darkfield):
    # correction by flat or dark field
    #calculate mean and std of the image to be able to normalize
    nbslices, height, width = Im.shape
    print (len(nbslices))
    meanSlice = np.mean(Im, axis=0)
    stdSlice = np.std(Im, axis=0)
    imCorrected = Im.copy()
    for slice in range(0, nbslices):
        imCorrected[slice,:,:]=(imCorrected[slice,:,:]-meanSlice)/stdSlice
        imCorrected[slice, :, :]= imCorrected[slice,:,:]-darkfield

    print ('-----------------------  normalization correction done ----------------------- ')
    return imCorrected







def normalizationMultipleDarkField (Im, darkfield):
    # correction by flat or dark field
    #calculate mean and std of the image to be able to normalize
    nbslices, height, width = Im.shape
    #meanSlice = np.mean(Im, axis=0)
    #meanSlice=gaussian_filter(meanSlice,sigma=50)
    #stdSlice = np.std(Im, axis=0)
    imCorrected=Im.copy()

    for slice in range(0, nbslices):
        #imCorrected[slice,:,:]=(Im[slice,:,:]-(meanSlice))#/stdSlice
        imCorrected[slice,:,:]=(imCorrected[slice,:,:]-darkfield[slice,:,:])#/(meanSlice -darkfield[slice,:,:])

    print ('-----------------------  normalization correction done ----------------------- ')
    return imCorrected




def normalization2D (Im, darkfield):
    # correction by flat or dark field
    #calculate mean and std of the image to be able to normalize
    meanSlice= np.mean(Im)
    stdSlice=np.std(Im)
    ImCorrected=(Im-meanSlice)/stdSlice
    ImCorrected=ImCorrected-darkfield
    print('-----------------------  normalization correction done ----------------------- ')
    return ImCorrected



if __name__ == "__main__":
    print('test corrections start')



    Ir = spytIO.openImage('ref1-1.edf')
    Is = spytIO.openImage('samp1-1.edf')
    DF= spytIO.openImage('DF1-1.edf')


    correctedRef=registerRefAndSample(Ir,Is,1000)
    correctedRef=normalization(correctedRef,DF)
    refCorrectedOutputFileName = '/Users/helene/PycharmProjects/spytlab/output/correctedRef.edf'
    outputEdf = edf.EdfFile(refCorrectedOutputFileName, access='wb+')
    outputEdf.WriteImage({}, correctedRef)


    correctedSample=normalization(Is,DF)
    sampleCorrectedOutputFileName = '/Users/helene/PycharmProjects/spytlab/output/correctedSample.edf'
    outputEdf = edf.EdfFile(sampleCorrectedOutputFileName, access='wb+')
    outputEdf.WriteImage({}, correctedSample)

