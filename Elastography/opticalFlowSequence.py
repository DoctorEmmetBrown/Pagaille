import glob
from utils import spytMkDir as mkdir
import spytIO
import OpticalFlow
import numpy as np


def createFolders(output):
    dxFolder = output + '/dx/'
    dyFolder = output + '/dy/'
    phiFolder = output + '/phi/'
    phi2Folder = output + '/phiKottler/'
    phi3Folder = output + '/phiLarkin/'
    gradientNormFolder=output + '/gradientNorm/'
    diffFolder = output + '/diff/'


    mkdir(dxFolder)
    mkdir(dyFolder)
    mkdir(phiFolder)
    mkdir(phi2Folder)
    mkdir(phi3Folder)
    mkdir(gradientNormFolder)
    mkdir(diffFolder)


def saveResults(result,outputFolder,numeroProj):
    dx = result['dx']
    dy = result['dy']
    phi = result['phi']
    phi2 = result['phi2']
    phi3 = result['phi3']
    gradientNorm = result['gradientNorm']

    textProj = '%4.4d' % numeroProj
    spytIO.saveEdf(dx.real, outputFolder + '/dx/dx_' + textProj + '.edf')
    spytIO.saveEdf(dy.real, outputFolder + '/dy/dy_' + textProj + '.edf')
    spytIO.saveEdf(phi.real, outputFolder  + '/phi/phi_' + textProj + '.edf')
    spytIO.saveEdf(phi2.real, outputFolder + '/phiKottler/phiKottler_' + textProj + '.edf')
    spytIO.saveEdf(phi3.real, outputFolder + '/phiLarkin/phiLarkin_' + textProj + '.edf')
    spytIO.saveEdf(gradientNorm.real, outputFolder  + '/gradientNorm/gradientNorm_' + textProj + '.edf')




if __name__ == "__main__":
    originalFolder='/Users/embrun/Downloads/Mantela_aurantiaca_1/cine1/'
    fileSequence=glob.glob(originalFolder+'/*.tif')
    fileSequence.sort()
    outputFolder='/Users/embrun/Downloads/Mantela_aurantiaca_1/Mantela_cine_OpticalFlow/'
    mkdir(outputFolder)
    createFolders(outputFolder)



    Ir=spytIO.openImage(fileSequence[0])
    Ir=np.asarray(Ir,np.float32)
    firstIm= spytIO.openImage(fileSequence[0])

    cpt=0
    for filename in fileSequence:
        print(filename)
        if cpt>0 :
            Is=spytIO.openImage(fileSequence[cpt])
            Is=np.asarray(Is,np.float32)
            result=OpticalFlow.processOneProjection(Is,Ir)
            saveResults(result,outputFolder,cpt)

            textProj = '%4.4d' % cpt
            spytIO.saveEdf(Is-Ir, outputFolder + '/diff/diff' + textProj + '.edf')
            Ir = Is

        cpt+=1


