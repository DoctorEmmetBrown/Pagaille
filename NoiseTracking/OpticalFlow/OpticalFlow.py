from InputOutput.pagailleIO import saveEdf,openImage,openSeq
import glob
import os
import sys
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2
from numpy.fft import fftfreq as fftfreq

from scipy.ndimage.filters import gaussian_filter
from math import pi as pi
from math import floor as floor

from NoiseTracking.OpticalFlow import frankoChellappa as fc
#import spytlabQT as qt
#import corrections

import numpy as np
from Tomography.FastTomo import fastTomoExperiment as esrfTomo


def derivativesByOpticalflow(intensityImage,derivative,alpha=0,sig_scale=0):

    Nx, Ny = derivative.shape #Image dimentions
    # fourier transfomm of the derivative and shift low frequencies to the centre
    ftdI = fftshift(fft2(derivative)) #Fourier transform of the derivative
    # calculate frequencies
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx) #frequency ranges of the images in fqcy space


    #building filters


    sigmaX = dqx / 1. * np.power(sig_scale,2) #0
    sigmaY = dqy / 1. * np.power(sig_scale,2)
    #sigmaX=sig_scale
    #sigmaY = sig_scale

    g = np.exp(-(((Qx)**2) / 2. / sigmaX + ((Qy)**2) / 2. / sigmaY))
    #g = np.exp(-(((np.power(Qx, 2)) / 2) / sigmaX + ((np.power(Qy, 2)) / 2) / sigmaY))
    beta = 1 - g;

    # fourier filters
    ftfiltX = (1j * Qx / ((Qx**2 + Qy**2 + alpha))*beta)
    ftfiltX[np.isnan(ftfiltX)] = 0

    ftfiltY = (1j* Qy/ ((Qx**2 + Qy**2 + alpha))*beta)
    ftfiltY[np.isnan(ftfiltY)] = 0

    # output calculation
    dImX = 1. / intensityImage * ifft2(ifftshift(ftfiltX * ftdI)) #Displacement field
    dImY = 1. / intensityImage * ifft2(ifftshift(ftfiltY * ftdI))

    return dImX.real,dImY.real




def kottler(dX,dY):
    print('kottler')
    i = complex(0, 1)
    Nx, Ny = dX.shape
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)

    polarAngle = np.arctan2(Qx, Qy)
    ftphi = fftshift(fft2(dX + i * dY))*np.exp(i*polarAngle)
    ftphi[np.isnan(ftphi)] = 0
    phi3 = ifft2(ifftshift(ftphi))
    return phi3.real



def LarkinAnissonSheppard(dx,dy,alpha =0 ,sigma=0):
    Nx, Ny = dx.shape
    i = complex(0, 1)
    G= dx + i*dy
    # fourier transfomm of the G function
    fourrierOfG = fftshift(fft2(G))


    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)

    ftfilt = 1 / (i * Qx - Qy)
    ftfilt[np.isnan(ftfilt)] = 0
    phi=ifft2(ifftshift(ftfilt*fourrierOfG))
    phi=phi.real
    return phi





def processOneProjection(Is,Ir):
    sigma = 0.9
    alpha = 0

    dI = (Is - Ir * (np.mean(gaussian_filter(Is,sigma=2.)) / np.mean(gaussian_filter(Ir,sigma=2.))))
    alpha=np.finfo(np.float32).eps
    dx, dy = derivativesByOpticalflow(Is, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    mintoAdd = min(np.amin(dx), np.amin(dy))
    mintoAdd=1
    dxPlus = dx + mintoAdd
    dyPlus = dy + mintoAdd

    gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm}


def processOneProjectionNormalizedSansWhiteField(Is,Ir,df):
    sigma = 0.9
    alpha = 0

    Isnormalized=(Is-df)/(Ir-df)
    Isnormalized[np.isnan(Isnormalized)]=0.000000001
    #saveEdf(Isnormalized,'output/Isnormalized.edf')
    Irnomralized= (Ir - df)/Ir

    dI = ((Isnormalized - Irnomralized)* (np.mean(gaussian_filter(Isnormalized,sigma=2.)) / np.mean(gaussian_filter(Irnomralized,sigma=2.))))

    #saveEdf(dI, 'output/dIdt.edf')
    alpha=np.finfo(np.float32).eps
    dx, dy = derivativesByOpticalflow(Irnomralized, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    mintoAdd = min(np.amin(dx), np.amin(dy))
    mintoAdd=5
    dxPlus = dx + mintoAdd
    dyPlus = dy + mintoAdd

    gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm}



def processOneProjectionNormalizedSansWhiteField2(Is,Ir,df):
    sigma = 0.9
    alpha = 0

    Isnormalized=(Is-df)
    spytIO.saveEdf(Isnormalized,'output/Isnormalized.edf')
    Irnomralized= (Ir - df)

    dI = ((Isnormalized - Irnomralized)* ((gaussian_filter(Isnormalized,sigma=2.)) / (gaussian_filter(Irnomralized,sigma=2.))))

    spytIO.saveEdf(dI, 'output/dIdt.edf')
    alpha=np.finfo(np.float32).eps
    dx, dy = derivativesByOpticalflow(Irnomralized, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    mintoAdd = min(np.amin(dx), np.amin(dy))
    mintoAdd=5
    dxPlus = dx + mintoAdd
    dyPlus = dy + mintoAdd

    gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm}





def processOneProjectionNormalized(Is,Ir,wf,df):
    sigma = 0.9
    alpha = 0

    Isnormalized=np.true_divide((Is-df),(wf-df))
    spytIO.saveEdf(Isnormalized, 'output/Isnormalized.edf')
    Irnomralized = np.true_divide((Ir - df), (wf - df))
    spytIO.saveEdf(Irnomralized,'output/Irnormalized.edf')

    dI = ((Isnormalized - Irnomralized) * (np.mean(gaussian_filter(Isnormalized,sigma=2.)) / np.mean(gaussian_filter(Irnomralized,sigma=2.))))
    alpha=np.finfo(np.float32).eps
    dx, dy = derivativesByOpticalflow(Irnomralized, dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)

    mintoAdd = min(np.amin(dx), np.amin(dy))
    mintoAdd=5
    dxPlus = dx + mintoAdd
    dyPlus = dy + mintoAdd

    gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm}



def processProjectionSetWithDarkFields(Is,Ir,dark):
    sigma = 1
    alpha = 0

    Is=corrections.normalizationMultipleDarkField(Is,dark)
    Ir=corrections.normalizationMultipleDarkField(Ir,dark)

    subImage=Is-Ir
    subImage=np.median(subImage,axis=0)

    dI = (subImage * (np.mean(Is) / np.mean(Ir)))
    dx, dy = derivativesByOpticalflow(np.mean(Ir,axis=0), dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)
    mintoAdd = min(np.amin(dx), np.amin(dy))
    dxPlus = dx + mintoAdd
    dyPlus = dy + mintoAdd

    gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm}


def processProjectionSetMedian(Is,Ir):
    sigma = 0.9
    alpha = 0

    subImageSet=Is-Ir
    print('processProjectionSetMedian')
    print(Is.shape)
    nbCoupleOfPoints=Is.shape[0]
    print(nbCoupleOfPoints)
    dxSet=np.zeros_like(Is,dtype=np.float32)
    dySet=np.zeros_like(Is, dtype=np.float32)
    phiSet=np.zeros_like(Is, dtype=np.float32)
    phi2Set=np.zeros_like(Is, dtype=np.float32)
    phi3Set=np.zeros_like(Is, dtype=np.float32)
    gradientNormSet=np.zeros_like(Is, dtype=np.float32)

    for i in range(0,nbCoupleOfPoints):
        subImage=subImageSet[i,:,:]
        dI = (subImage * (np.mean(gaussian_filter(Is[i,:,:],sigma=2.)) / np.mean(gaussian_filter(Ir[i,:,:],sigma=2.))))
        dx, dy = derivativesByOpticalflow(np.mean(Ir, axis=0), dI, alpha=alpha, sig_scale=sigma)
        dxSet[i,:,:]=dx
        dySet[i, :, :]=dy
        phi = fc.frankotchellappa(dx, dy, False)
        phi3 = kottler(dx, dy)
        phi2 = LarkinAnissonSheppard(dx, dy)
        mintoAdd = min(np.amin(dx), np.amin(dy))
        dxPlus = dx + mintoAdd
        dyPlus = dy + mintoAdd
        gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)
        phiSet[i,:,:]=phi.real
        phi2Set[i, :, :] = phi2.real
        phi3Set[i, :, :] = phi3.real
        gradientNormSet[i, :, :] = gradientNorm

    dx=np.median(dxSet,axis=0)
    dy = np.median(dySet, axis=0)
    phi= np.median(phiSet, axis=0)
    phi2 = np.median(phi2Set, axis=0)
    phi3 = np.median(phi3Set, axis=0)
    gradientNorm = np.median(gradientNormSet, axis=0)
    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm }

def processProjectionSetMedian2(Is,Ir):
    sigma = 0.9
    alpha = 0

    subImageSet=Is-Ir
    print('processProjectionSetMedian')
    print(Is.shape)
    nbCoupleOfPoints=Is.shape[0]
    print(nbCoupleOfPoints)
    dxSet=np.zeros_like(Is,dtype=np.float32)
    dySet=np.zeros_like(Is, dtype=np.float32)
    phiSet=np.zeros_like(Is, dtype=np.float32)
    phi2Set=np.zeros_like(Is, dtype=np.float32)
    phi3Set=np.zeros_like(Is, dtype=np.float32)
    gradientNormSet=np.zeros_like(Is, dtype=np.float32)

    for i in range(0,nbCoupleOfPoints):
        subImage=subImageSet[i,:,:]
        dI = (subImage * (np.mean(gaussian_filter(Is[i,:,:],sigma=2.)) / np.mean(gaussian_filter(Ir[i,:,:],sigma=2.))))
        dx, dy = derivativesByOpticalflow(np.mean(Ir, axis=0), dI, alpha=alpha, sig_scale=sigma)
        dxSet[i,:,:]=dx
        dySet[i, :, :]=dy

    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)
    mintoAdd = min(np.amin(dx), np.amin(dy))
    dxPlus = dx + mintoAdd
    dyPlus = dy + mintoAdd

    gradientNorm = np.sqrt(dxPlus ** 2 + dyPlus ** 2)

    dx=np.median(dxSet,axis=0)
    dy = np.median(dySet, axis=0)




    gradientNorm = np.median(gradientNormSet, axis=0)
    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm }






def processProjectionSet(Is,Ir):
    sigma = 0.9
    alpha = 0

    subImage=Is-Ir
    subImage=np.median(subImage,axis=0)

    dI = (subImage * (np.mean(gaussian_filter(Is,sigma=2)) / np.mean(gaussian_filter(Ir,sigma=2))))
    dx, dy = derivativesByOpticalflow(np.mean(Ir,axis=0), dI, alpha=alpha, sig_scale=sigma)
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)
    mintoAdd = min(np.amin(dx), np.amin(dy))
    dxPlus=dx+mintoAdd
    dyPlus=dy+mintoAdd

    gradientNorm = np.sqrt(dxPlus**2 + dyPlus**2)

    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3,'gradientNorm':gradientNorm }


if __name__ == "__main__":    #Ir = spytIO.openImage('/Users/embrun/Codes/specklematching/Experiments/Test_id17_fev2018_fils_2D/ref/Ref0073.edf')
    #Is = spytIO.openImage('/Users/embrun/Codes/specklematching/Experiments/Test_id17_fev2018_fils_2D/sample/im0001.edf')
    #wf=np.zeros_like(Ir)
    #df = np.zeros_like(Ir)


    #Ir = spytIO.openImage('/Volumes/VISITOR/md1176/id17/HA500_Speckle_Dent_34kev_6um_bis__000__003_/refHST3000.edf')
    #Is = spytIO.openImage('/Volumes/VISITOR/md1176/id17/HA500_Speckle_Dent_34kev_6um_bis__000__003_/HA500_Speckle_Dent_34kev_6um_bis__000__003_0002.edf')
    #df = spytIO.openImage('/Volumes/VISITOR/md1176/id17/HA500_Speckle_Dent_34kev_6um_bis__000__003_/dark.edf')
    #wf=gaussian_filter(Ir,6.)
    #wf = spytIO.openImage('/Volumes/VISITOR/md1176/id17/HA500_Speckle_Dent_34kev_6um_bis__000__Propag_/refHST3000.edf')
    #wf = gaussian_filter(Ir, 3.)

    Ir = openImage('Z:/speckle2/TestLadafMarch2019TestImpression/Tif/PremiersTests/BalleRef/ESRF.XA._.5.1.2019.03.20.15.14.58.671875.4038525.tif')
    Is = openImage('Z:/speckle2/TestLadafMarch2019TestImpression/Tif/PremiersTests/Balle/ESRF.XA._.4.1.2019.03.20.15.14.58.671875.4038488.tif')
    Ir = Ir[56:56 + 570, 414:414 + 600]
    Is = Is[56:56 + 570, 414:414 + 600]


    #Ir = Ir[1061:1061+ 440, 915:915+ 700]
    #Is = Is[1061:1061+ 440, 915:915+ 700]

    #Ir=registerRefAndSample(Ir, Is, 1000)

    #df = spytIO.openImage('/Volumes/VISITOR/md1176/id17/HA500_Speckle_Dent_34kev_6um_bis__000__003_/dark.edf')
    #wf=gaussian_filter(Ir,6.)
    #wf = spytIO.openImage('/Volumes/VISITOR/md1176/id17/HA500_Speckle_Dent_34kev_6um_bis__000__Propag_/refHST3000.edf')
    #wf = gaussian_filter(Ir, 3.)



#    Ir = spytIO.openImage('/Volumes/ID17/speckle2/md1125/id17/Foie/NASH/NASSH_16PH02530_SpeckleCu45_11um_al4_cu045__003_/refForHST2000.edf')
#    Is = spytIO.openImage('/Volumes/ID17/speckle2/md1125/id17/Foie/NASH/NASSH_16PH02530_SpeckleCu45_11um_al4_cu045__003_/NASSH_16PH02530_SpeckleCu45_11um_al4_cu045__003_0008.edf')
#    df = spytIO.openImage('/Volumes/ID17/speckle2/md1125/id17/Foie/NASH/NASSH_16PH02530_SpeckleCu45_11um_al4_cu045__003_/darkForHST0000.edf')


    #Ir = spytIO.openImage('/Volumes/ID17/speckle2/md1125/id17/Tresse/Tresse_Speckle_Cu45_11um_PINK_1800prj_6pts__001_/ref0006_1800.edf')
    #Is = spytIO.openImage('/Volumes/ID17/speckle2/md1125/id17/Tresse/Tresse_Speckle_Cu45_11um_PINK_1800prj_6pts__001_/Tresse_Speckle_Cu45_11um_PINK_1800prj_6pts__001_0019.edf')
    #df = spytIO.openImage('/Volumes/ID17/speckle2/md1125/id17/Tresse/Tresse_Speckle_Cu45_11um_PINK_1800prj_6pts__001_/darkend0000.edf')
    #df = df/10.



    #Is = spytIO.openImage('/Volumes/ID17/speckle/DonneeCArmLadaf/Tiff/PolentaPoulet.tif')
    #Ir = spytIO.openImage('/Volumes/ID17/speckle/DonneeCArmLadaf/Tiff/PolentaReference.tif')

    #Ir = spytIO.openImage('/Volumes/ID17/inhouse4/MITTONE/md1127/AK17-5639/SPECKLES/AK175639_speckles_BIN4_1ms_360p_4cm_polenta_/AK175639_speckles__BIN4_1ms_360p_4cm__004_1/AK175639_speckles__BIN4_1ms_360p_4cm__004_10018.edf')
    #Is = spytIO.openImage('/Volumes/ID17/inhouse4/MITTONE/md1127/AK17-5639/SPECKLES/AK175639_speckles_BIN4_1ms_360p_4cm_polenta_/AK175639_speckles__BIN4_1ms_360p_4cm__004_1/refForHST0000.edf')
    #df = spytIO.openImage('/Volumes/ID17/inhouse4/MITTONE/md1127/AK17-5639/SPECKLES/AK175639_speckles_BIN4_1ms_360p_4cm_polenta_/AK175639_speckles__BIN4_1ms_360p_4cm__004_1/darkend0000.edf')


    #Ir=Ir[264:264+460,65:65+839]
    #Is=Is[264:264+460,65:65+839]
 #   df = np.random.rand(Ir.shape)
    #Is=Is+1
    #Ir=Ir+1

    #print(Ir.shape)

    #header=spytIO.getHeader('/Users/embrun/Codes/specklematching/Experiments/MoucheSimapAout2017/ref/ref_40kV_5.2um_30s_12cm_53cm_speck06.tif')

    #wf = wf[86:86 + 294, 82:82 + 2410]
    #df = df[86:86 + 294, 82:82 + 2410]
    nrow,ncol=Ir.shape
    print(nrow)
    print(ncol)
    padCol=nrow
    padRow=ncol

    Ir=np.pad(Ir,((padRow,padRow),(padCol,padCol)),'reflect')
    Is=np.pad(Is,((padRow,padRow),(padCol,padCol)), 'reflect')
    #df = np.pad(df, ((padRow, padRow), (padCol, padCol)), 'reflect')
    #wf = np.pad(wf, ((padRow, padRow), (padCol, padCol)), 'reflect')
    #Ir = Ir[86:86+294,82:82+ 2410]
    #Is=Is[86:86+294,82:82+ 2410]


    Ir=np.asarray(Ir,dtype=np.float32)
    Is=np.asarray(Is, dtype=np.float32)
    #df = np.asarray(df, dtype=np.float32)


    #result = processOneProjection(Is, Ir)
    #result = processOneProjectionNormalized(Is, Ir,wf,df)
    #result = processOneProjectionNormalizedSansWhiteField(Is, Ir,df)

    #IrNames=glob.glob('/Users/embrun/Codes/specklematching/Experiments/projection_WireAndSpheres/ref//*.tiff')
    #IsNames= glob.glob('/Users/embrun/Codes/specklematching/Experiments/projection_WireAndSpheres/sample/*.tiff')
    #Ir=spytIO.openSeq(IrNames)
    #Is= spytIO.openSeq(IsNames)
    #result = processProjectionSetMedian2(Is, Ir)

    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    gradientNorm = result['gradientNorm']
    Ir=np.asarray(Ir,dtype=np.float32)
    Is=np.asarray(Is,dtype=np.float32)

    Is=Is[padRow:padRow+nrow,padCol:padCol+ncol]
    Ir=Ir[padRow:padRow+nrow,padCol:padCol+ncol]
    dx = dx[padRow:padRow+nrow,padCol:padCol+ncol]
    dy = dy[padRow:padRow+nrow,padCol:padCol+ncol]
    phi = phi[padRow:padRow+nrow,padCol:padCol+ncol]
    phi2 = phi2[padRow:padRow+nrow,padCol:padCol+ncol]
    phi3 = phi3[padRow:padRow+nrow,padCol:padCol+ncol]
    gradientNorm= gradientNorm[padRow:padRow+nrow,padCol:padCol+ncol]
    tmp=Ir-Is

    outputFolder='C:/Users/quenot/Desktop/outputOpticalFlow/Ball'
    saveEdf(dx, outputFolder+'/dx.edf')
    saveEdf(dy.real, outputFolder+'/dy.edf')
    saveEdf(phi.real, outputFolder+'/phi.edf')
    saveEdf(phi2.real, outputFolder+'/phiLarkinson.edf')
    saveEdf(phi3.real, outputFolder+'/phiKottler.edf')
    saveEdf(gradientNorm, outputFolder+'/gradientNorm.edf')
    saveEdf(tmp, outputFolder+'/soustraction.edf')
    #spytIO.saveEdf(Is, 'output/sample.edf')
    #spytIO.saveTiff16bit(gradientNorm, 'output/gradientNorm.tif',header=header)
