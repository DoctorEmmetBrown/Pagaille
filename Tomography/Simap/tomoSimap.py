import numpy as np
import corrections as corr
import fastTomoExperiment as esrfTomo
import OpticalFlow
import spytlabQT as qt
import glob
import os
import spytIO
import opticalThread
from  utils import spytMkDir as mkdir
import Averager


def opticalFlowSolverOneThread(inputDictionary, outputFolder):
    projectionFiles = inputDictionary['projections']
    referenceFilename = inputDictionary['references'][0]
    Ir = spytIO.openImage(referenceFilename)
    i = 0
    for projectionFileName in projectionFiles:
        print(projectionFileName)
        Is = spytIO.openImage(projectionFileName)
        result = OpticalFlow.processOneProjection(Is, Ir)
        saveResult(result, i, outputFolder)
        i += 1


def saveResult(res, projNumber, outputFolder):
    dx = res['dx']
    dy = res['dy']
    phi = res['phi']
    phi2 = res['phi2']
    phi3 = res['phi3']
    txtProj = '%4.4d' % projNumber
    spytIO.saveEdf(dx, outputFolder + '/dx/dx_' + txtProj + '.edf')
    spytIO.saveEdf(dy.real, outputFolder + '/dy/dy_' + txtProj + '.edf')
    spytIO.saveEdf(phi.real, outputFolder + '/phi/phi_' + txtProj + '.edf')
    spytIO.saveEdf(phi2.real, outputFolder + '/phiLarkin/phiLarkin_' + txtProj + '.edf')
    spytIO.saveEdf(phi3.real, outputFolder + '/phiKottler/phiKottler' + txtProj + '.edf')


def parseTomoFolderAndCreateRefFiles(folderpath):
    scanName = os.path.basename(folderpath)

    projectionsFileNames=glob.glob(folderpath+'tomo_point1*')
    projectionsFileNames.sort()
    nbProjection=len(projectionsFileNames)
    txtNbProj='%4.4d'%nbProjection

    referenceApresFileNames=glob.glob(folderpath+'/ref_apres*')
    referenceApresFileNames.sort()


    referenceAvantFileNames = glob.glob(folderpath + '/ref_apres*')

    referenceAvantFileNames.sort()


    referenceFileNames = glob.glob(folderpath + '/*refForHST*')
    if referenceFileNames == []:
        outputFileFFNameend = folderpath + '/refForHST' + txtNbProj + '.edf'

        Averager.Averager(referenceApresFileNames, outputFileFFNameend, option=0)
        outputFileFFNameBeg = folderpath + '/refForHST0000.edf'
        Averager.Averager(referenceAvantFileNames, outputFileFFNameBeg, option=0)
        referenceFileNames = glob.glob(folderpath + '/*refForHST*')

    darkFieldFilename = 'bidon.edf'
    print(darkFieldFilename)
    referenceFileNames.sort()
    print(referenceFileNames)
    ddict = {}
    ddict['projections'] = projectionsFileNames
    ddict['references'] = referenceFileNames
    ddict['darkField'] = darkFieldFilename
    return ddict


def createOutput(outputfolder):
    mkdir(outputfolder)
    dxFolder = outputfolder + '/dx/'
    dyFolder = outputfolder + '/dy/'
    phiFolder = outputfolder + '/phi/'
    phi2Folder = outputfolder + '/phiLarkin/'
    phi3Folder = outputfolder + '/phiKottler/'
    mkdir(dxFolder)
    mkdir(dyFolder)
    mkdir(phiFolder)
    mkdir(phi2Folder)
    mkdir(phi3Folder)


def opticalFlowSolverMultipleThread(inputFolder, outputFolder, nbThread=4):
    ddict = parseTomoFolderAndCreateRefFiles(inputFolder)
    numberOfProjections = len(ddict['projections'])
    listofThreads = []
    nbProjByThread = int(numberOfProjections / nbThread)
    print('nbProjByThread' + str(nbProjByThread))
    for i in range(nbThread):
        if i == nbThread - 1:
            listOfProjections = (np.arange(i * nbProjByThread, numberOfProjections))
        else:
            listOfProjections = (np.arange(i * nbProjByThread, (i + 1) * nbProjByThread))

        myThread = opticalThread.OpticalFlowSolverThread(ddict, listOfProjections, outputFolder)
        listofThreads.append(myThread)

    for i in range(nbThread):
        listofThreads[i].start()

    for i in range(nbThread):
        listofThreads[i].join()


if __name__ == "__main__":
    inputFolder = '/Volumes/KINGSTON/proj_fin_cigale/'
    outputFolder = '/Users/helene/Desktop/OutputCigalePhi'
    createOutput(outputFolder)
    ddict=parseTomoFolderAndCreateRefFiles(inputFolder)
    opticalFlowSolverMultipleThread(inputFolder,outputFolder,nbThread=4)


