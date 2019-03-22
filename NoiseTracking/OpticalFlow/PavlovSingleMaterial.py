from InputOutput.pagailleIO import openImage,saveEdf
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2
from math import pi as pi
import numpy as np
import os
import glob
from NoiseTracking.OpticalFlow import pavlovThread



def tie_hom_KMP2(Im, Energy, sourceToSampleDist, sampleToDetDist, pix_size, delta, beta):
    waveNumber=kev2wn(Energy)

    factor = magnificationFactor(sourceToSampleDist, sampleToDetDist)
    print('magnification factor ' + str(factor))

    factor=1
    Im=Im*(factor**2)
    nrow,ncol=Im.shape
    padRow=1
    padCol=1
    Im = np.pad(Im, ((padRow, padRow), (padCol, padCol)), 'reflect')
    ImFft=fftshift(fft2(Im))

    Nx, Ny = ImFft.shape
    # calculate frequencies
    kx, ky = np.meshgrid(np.arange(0, Nx), np.arange(0, Ny))
    kx = (kx - (Nx / 2))
    ky = (ky - (Ny / 2))
    kx = 2 * pi * kx / (Nx * pix_size)
    ky = 2 * pi * ky / (Ny * pix_size)

    k = np.sqrt(kx**2+ ky** 2)
    k_sqr = np.transpose(k) **2

    ffilt = k_sqr * sampleToDetDist * factor * delta + 2 * waveNumber * beta
    tmp = ImFft/ ffilt
    tmp =ifft2(ifftshift(tmp))

    ImThickness = tmp.real
    ImThickness = ImThickness[padRow:padRow + nrow, padCol:padCol + ncol]

    return ImThickness


def tie_hom_KMP2Last (img_in, E, R1, R2, pix_size, RI_delta, RI_beta, bg_val, scale):
    """ Phase retrieval using the Transport of Intensity Equation for homogeneous samples in the case of speckles
     Parameters
    ----------
    arg1 : img_in
       Image to be phase retrieved (must be reference corrected)
    arg2 : str
        E = mean x-ray energy (keV)
    arg3: float
        R1 = source to sample distance (m)
    arg4: float
        R2 = sample to detector distance (m)
    arg5: float
        pix_size = size of the image pixels (m)

    arg6: float
        RI_delta = refractive index

    arg7: float
        RI_beta = absorption par    t of the refactive index

    arg8: float
        bg_val = background intensity in img_in (e.g. 0, 1, 100...)

    arg9: float
        scale = parameter in the low-frequency filter (e.g., 2) (added by KMP)

    Returns
    -------
    float
        img_thickness = image of sample thickness (m)

    """

    k_in = kev2wn(E)
    magnificationFactor=(R1+R2)/R1
    padCol = 600
    padRow= 600

    width,height=img_in.shape
    img_in = np.pad(img_in, ((padRow, padRow), (padCol, padCol)), 'reflect')

    #img_obj_plane = np.zeros_like(img_in)
    #img_thickness = img_obj_plane
    img_infft = fftshift(fft2(img_in))
    Nx, Ny = img_infft.shape
    kx, ky = np.meshgrid(np.arange(0, Nx), np.arange(0, Ny))
    kx = (kx - (Nx / 2))
    ky = (ky - (Ny / 2))
    kx = 2 * pi * kx / (Nx * pix_size)
    ky = 2 * pi * ky / (Ny * pix_size)

    k = np.sqrt(kx ** 2 + ky ** 2)
    k_sqr = np.transpose(k) ** 2

    sigma_x = ((2 * pi / (Nx * pix_size)) * scale)** 2
    sigma_y = ((2 * pi / (Ny * pix_size)) * scale) ** 2
    low_freq = (1. - np.exp(-(kx**2 / (2. * sigma_x) + ky**2 / (2. * sigma_y))))
    low_freq1 = np.transpose(low_freq)

    #dbratio = RI_delta / RI_beta

    ffilt = k_sqr * R2 * RI_delta / magnificationFactor + 2 * k_in * RI_beta

    img_infft = img_infft * low_freq1

    tmp = img_infft / ffilt
    tmp = (ifft2(ifftshift(tmp)))
    img_thickness = np.real(tmp)
    img_thickness=img_thickness[padRow:padRow+width,padCol:padCol+height]

    return img_thickness


def tie_hom_KMP2_normalizexp(Im, Energy, sourceToSampleDist, sampleToDetDist, pix_size, delta, beta,sig_scale=0):
    waveNumber=kev2wn(Energy)

    factor = magnificationFactor(sourceToSampleDist, sampleToDetDist)
    print('magnification factor ' + str(factor))

    factor=1
    Im=Im*(factor**2)
    nrow,ncol=Im.shape
    padRow=1000
    padCol=1000
    Im = np.pad(Im, ((padRow, padRow), (padCol, padCol)), 'reflect')
    ImFft=fftshift(fft2(Im))

    Nx, Ny = ImFft.shape
    # calculate frequencies
    kx, ky = np.meshgrid(np.arange(0, Nx), np.arange(0, Ny))
    kx = (kx - (Nx / 2))
    ky = (ky - (Ny / 2))
    kx = 2 * pi * kx / (Nx * pix_size)
    ky = 2 * pi * ky / (Ny * pix_size)

    k = np.sqrt(kx**2+ ky** 2)
    k_sqr = np.transpose(k) **2


  #  sigmaX = Nx / 1. * np.power(sig_scale,2)
   # sigmaY = Ny / 1. * np.power(sig_scale,2)
    #g = np.exp(-(((kx)**2) / 2. / sigmaX + ((ky)**2) / 2. / sigmaY))
    #g = np.exp(-(((np.power(Qx, 2)) / 2) / sigmaX + ((np.power(Qy, 2)) / 2) / sigmaY))
    #beta = 1 - g
#    beta=np.transpose(beta)

    ffilt = k_sqr * sampleToDetDist * factor * delta + 2 * waveNumber * beta
    tmp = ImFft/ ffilt
    tmp =ifft2(ifftshift(tmp))

    ImThickness = tmp.real
    ImThickness = ImThickness[padRow:padRow + nrow, padCol:padCol + ncol]

    return ImThickness

def kev2wn(kev):
    #Converts energy in keV to corresponding wave number
    wn = kev*1000*2*pi*1.60221773e-19/(6.6260755e-34 * 2.99792458e8)
    print(wn)
    return wn

def magnificationFactor(sourceToSampleDist,sampleToDetDist):
    #Calculate geometrical magnification
    factor=(sourceToSampleDist+sampleToDetDist)/sourceToSampleDist
    return factor




def tomoPavlov(samplefilenames,referencename,outputFolder):
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

    Ir = openImage(referencename)
    Ir = np.asarray(Ir, np.float32)

    for filename in samplefilenames:
        print(filename)
        Is = openImage(filename)
        Is = np.asarray(Is, np.float32)

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

        img_out = tie_hom_KMP2Last(img_in, E, R1, R2, pix_size, delta, beta * low_frequency_filter, bg_val, scale)
        img_out = img_out * 1e6
        saveEdf(img_out,outputFolder+'/'+os.path.basename(filename))





def tomoPavlov(samplefilenames,referencename,outputFolder):
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

    Ir = openImage(referencename)
    Ir = np.asarray(Ir, np.float32)

    for filename in samplefilenames:
        print(filename)
        Is = openImage(filename)
        Is = np.asarray(Is, np.float32)

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

        img_out = tie_hom_KMP2Last(img_in, E, R1, R2, pix_size, delta, beta * low_frequency_filter, bg_val, scale)
        img_out = img_out * 1e6
        saveEdf(img_out,outputFolder+'/'+os.path.basename(filename))



def processOneImage(sampleImageName,referenceImageName,darkFileName=None,whiteFieldFileName=None):
    # Constants
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


    Ir = openImage(referenceImageName)
    Is = openImage(sampleImageName)
    Ir = np.asarray(Ir, np.float32)
    Is = np.asarray(Is, np.float32)

    if darkFileName is not None:
        df=openImage(darkFileName)
        Is = (Is - df) / (Ir - df)
        Is[np.isnan(Is)] = 0.0000000001
        Ir= (Ir - df) / Ir
        Is[np.isnan(Is)] = 0.0000000001


    Image_old = np.true_divide(Is, Ir)
    Image_old = 1 - Image_old

    # New smaller array containing image of fibres only

    # Copying a part of the larger array into the smaller array
    # Image_new = Image_old[2:900,2:2000]
    Image_new = Image_old

    # Calculation of the average value of the image
    average_image = np.mean(Image_new)
    # Correction on the average image. Now the average of the new array is ~0
    Image_new = Image_new - average_image
    saveEdf(Image_new, 'ImageNew.edf')

    # Using a modified TIE-HOM retrieval program from
    # https://github.com/RSBradley/TomoTools/blob/master/addins/phase%20retrieval/tie_hom.m

    img_in = Image_new

    E = 52
    low_frequency_filter = 1.
    scale = 20.

    # img_out = tie_hom_KMP2(img_in, E, R1, R2, pix_size, delta, beta*low_frequency_filter)
    bg_val = 0

    img_out = tie_hom_KMP2Last(img_in, E, R1, R2, pix_size, delta, beta * low_frequency_filter, bg_val, scale)
    # img_out = tie_hom_KMP2_normalizexp(img_in, E, R1, R2, pix_size, delta, beta,sig_scale=0.9)
    img_out = img_out * 1e6
    # img_out[img_out<0]=0.000000000001

    #saveEdf(img_out, '/Volumes/ID17/broncho/IHR_April2018/PavlovHA800Patte21_speckle01/HA800_Patte21_3um_Gap90_75_Speckle01_'+numSlice)

    return img_out



def tomoPavlovMultiThreaded(samplefilenames,referencename,outputFolder,darkName,nbThread=4):
    numberOfProjections=len(samplefilenames)
    listofThreads = []
    nbProjByThread = int(numberOfProjections / nbThread)
    print('nbProjByThread' + str(nbProjByThread))
    for i in range(nbThread):
        if i == nbThread - 1:
            listOfProjections = (np.arange(i * nbProjByThread, numberOfProjections))
        else:
            listOfProjections = (np.arange(i * nbProjByThread, (i + 1) * nbProjByThread))

        print(darkName)
        myThread = pavlovThread.PavlovOpticalFlowSolverThread(listOfProjections, samplefilenames,referencename,darkName,outputFolder)

        listofThreads.append(myThread)

    for i in range(nbThread):
        listofThreads[i].start()

    for i in range(nbThread):
        listofThreads[i].join()


def testOneImage():
    # Constants
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

    numSlice='0000.edf'

    # Reading data
    Ir = openImage('/VOLUMES/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/refHST0000.edf')
    Is = openImage('/VOLUMES/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/HA800_Patte21_3um_Gap90_75_Speckle01_'+numSlice)
    df=openImage('/VOLUMES/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/dark.edf')




    Ir = np.asarray(Ir, np.float32)
    Is = np.asarray(Is, np.float32)
    Is = (Is - df) / (Ir - df)
    Is[np.isnan(Is)] = 0.0000000001
    Ir= (Ir - df) / Ir
    Is[np.isnan(Is)] = 0.0000000001


    Image_old = np.true_divide(Is, Ir)
    Image_old = 1 - Image_old

    # New smaller array containing image of fibres only

    # Copying a part of the larger array into the smaller array
    # Image_new = Image_old[2:900,2:2000]
    Image_new = Image_old

    # Calculation of the average value of the image
    average_image = np.mean(Image_new)
    # Correction on the average image. Now the average of the new array is ~0
    Image_new = Image_new - average_image
    saveEdf(Image_new, 'ImageNew.edf')

    # Using a modified TIE-HOM retrieval program from
    # https://github.com/RSBradley/TomoTools/blob/master/addins/phase%20retrieval/tie_hom.m

    img_in = Image_new

    E = 52
    low_frequency_filter = 1.
    scale = 20.

    # img_out = tie_hom_KMP2(img_in, E, R1, R2, pix_size, delta, beta*low_frequency_filter)
    bg_val = 0

    img_out = tie_hom_KMP2Last(img_in, E, R1, R2, pix_size, delta, beta * low_frequency_filter, bg_val, scale)
    # img_out = tie_hom_KMP2_normalizexp(img_in, E, R1, R2, pix_size, delta, beta,sig_scale=0.9)
    img_out = img_out * 1e6
    # img_out[img_out<0]=0.000000000001

    #saveEdf(img_out, '/Volumes/ID17/broncho/IHR_April2018/PavlovHA800Patte21_speckle01/HA800_Patte21_3um_Gap90_75_Speckle01_'+numSlice)
    saveEdf(img_out, 'thickness.edf')




if __name__ == "__main__":
    #testOneImage()
    irName='/VOLUMES/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/refHST0000.edf'
    isName='/VOLUMES/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/HA800_Patte21_3um_Gap90_75_Speckle01_0000.edf'
    darkName='/VOLUMES/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/dark.edf'


    result=processOneImage(isName,irName,darkName)
    saveEdf(result, 'thickness.edf')
    outputtomofolder='/Volumes/ID17/broncho/IHR_April2018/PavlovHA800Patte21_speckle01/'
    inputFilenames=glob.glob('/Volumes/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/HA800*.edf')
    inputFilenames.sort()
    referenceFilename='/Volumes/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/refHST0000.edf'
    darkFileName= '/Volumes/ID17/broncho/IHR_April2018/HA800_Patte21_3um_Gap90_75_Speckle01_/dark.edf'
    tomoPavlovMultiThreaded(inputFilenames,referenceFilename,outputtomofolder,darkFileName)




