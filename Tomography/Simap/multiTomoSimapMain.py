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



def processOneProjection(listOfDictionnaries,projectionNumber):
    print('------------------------------------------------------')
    print('processOneProjection')
    projectionFiles=[]
    referencesFiles=[]
    darkFieldFiles=[]
    for dict in listOfDictionnaries :
        projectionFiles.append(dict['projections'][projectionNumber])
        referencesFiles.append(dict['references'][0])
        darkFieldFiles.append(dict['darkField'])

    print(projectionFiles)
    Is=spytIO.openSeq(projectionFiles)
    print(Is.shape)
    Ir=spytIO.openSeq(referencesFiles)
    print(Ir.shape)
    df=spytIO.openSeq(darkFieldFiles)
    print(df.shape)
    #Is,Ir=corr.registerImagesBetweenThemselves(Is,Ir)

    toReturn=OpticalFlow.processProjectionSetWithDarkFields(Is,Ir,df)
    return toReturn



def processAllFolders(listOfFolders,outputFolder):
    dxFolder = outputFolder + '/dx/'
    dyFolder = outputFolder + '/dy/'
    phiFolder = outputFolder + '/phi/'
    phi2Folder = outputFolder + '/phi2/'
    phi3Folder = outputFolder + '/phi3/'
    mkdir(dxFolder)
    mkdir(dyFolder)
    mkdir(phiFolder)
    mkdir(phi2Folder)
    mkdir(phi3Folder)


    listOfDictionaries=[]
    for folder in listOfFolders:
        ddict=parseTomoFolderAndCreateRefFiles(folder)
        listOfDictionaries.append(ddict)

    numberOfProjections=len(listOfDictionaries[0]['projections'])
    for projectionNumber in range (0,numberOfProjections):
        projectionNumber=1200
        result=processOneProjection(listOfDictionaries,projectionNumber)
        textProj='%4.4d'%projectionNumber

        dx = result['dx']
        dy = result['dy']
        phi = result['phi']
        phi2 = result['phi2']
        phi3 = result['phi3']
        spytIO.saveEdf(dx, dxFolder + '/dx'+textProj+'.edf')
        spytIO.saveEdf(dy.real, dyFolder + '/dy'+textProj+'.edf')
        spytIO.saveEdf(phi.real, phiFolder + '/phi_'+textProj+'.edf')
        spytIO.saveEdf(phi2.real, phi2Folder + '/phiKottler_'+textProj+'.edf')
        spytIO.saveEdf(phi3.real, phi3Folder + '/phiLarkin_'+textProj+'.edf')






def processAllFoldersThreaded(listOfFolders,outputFolder,nbThread=4):
    dxFolder = outputFolder + '/dx/'
    dyFolder = outputFolder + '/dy/'
    phiFolder = outputFolder + '/phi/'
    phi2Folder = outputFolder + '/phiKottler/'
    phi3Folder = outputFolder + '/phiLarkin/'
    mkdir(dxFolder)
    mkdir(dyFolder)
    mkdir(phiFolder)
    mkdir(phi2Folder)
    mkdir(phi3Folder)


    listOfDictionaries=[]
    for folder in listOfFolders:
        ddict=parseTomoFolderAndCreateRefFiles(folder)
        listOfDictionaries.append(ddict)

    numberOfProjections=len(listOfDictionaries[0]['projections'])

    listofThreads=[]
    nbProjByThread=int(numberOfProjections/nbThread)
    print('nbProjByThread'+str(nbProjByThread))
    for i in range(nbThread):
            listOfProjections = (np.arange(i * nbProjByThread,numberOfProjections))
        else:
            listOfProjections = (np.arange(i*nbProjByThread,(i+1)*nbProjByThread))

        myThread=opticalThread.multiTomoOpticalFlowSolver(listOfDictionaries,listOfProjections,outputFolder)
        listofThreads.append(myThread)

    for i in range(nbThread):
        listofThreads[i].start()

    for i in range(nbThread):
        listofThreads[i].join()





def parseTomoFolderAndCreateRefFiles(folderpath):
    scanName = os.path.basename(folderpath)

    projectionsFileNames=glob.glob(folderpath+'/mouche_speckle_*')
    projectionsFileNames.sort()
    print projectionsFileNames
    nbProjection=len(projectionsFileNames)
    txtNbProj='%4.4d'%nbProjection

    referenceApresFileNames=glob.glob(folderpath+'/ref_apres*')
    referenceApresFileNames.sort()


    referenceAvantFileNames = glob.glob(folderpath + '/ref_avant*')

    referenceAvantFileNames.sort()

    referenceFileNames = glob.glob(folderpath + '/*refForHST*')
    if referenceFileNames == [] :
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




if __name__ == "__main__":
    inputFolder='/Volumes/ID17/rsrm/Tavelures/SIMAP_112017/2017_11_08_TestSpeckle/tomo_mouche_speckle/'
    outputFolder = '/Volumes/ID17/rsrm/Tavelures/SIMAP_112017/2017_11_08_TestSpeckle/tomo_mouche_speckle/OpticalFlowMultiplePoints/'
    mkdir(outputFolder)
    tomoFolders=glob.glob(inputFolder+'renamed_tomo_mouche*')
    tomoFolders.sort()
    print(tomoFolders)
    processAllFoldersThreaded(tomoFolders,outputFolder,nbThread=4)

    #result=parseTomoFolderAndCreateRefFiles(tomoFolders[0])
    #print result


