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
        spytIO.saveEdf(phi.real, phiFolder + '/phi'+textProj+'.edf')
        spytIO.saveEdf(phi2.real, phi2Folder + '/phi2'+textProj+'.edf')
        spytIO.saveEdf(phi3.real, phi3Folder + '/phi3'+textProj+'.edf')






def processAllFoldersThreaded(listOfFolders,outputFolder,nbThread=4):
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

    listofThreads=[]
    nbProjByThread=int(numberOfProjections/nbThread)
    print('nbProjByThread'+str(nbProjByThread))
    for i in range(nbThread):
        if i == nbThread-1:
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
    scanName=os.path.basename(folderpath)
    parametersScanFilename=folderpath+'/'+scanName+'.xml'
    print(parametersScanFilename)
    tomoExperiment=esrfTomo.FastTomoExperiment(parametersScanFilename)
    print('numberFlatField: ')
    print(tomoExperiment.numberFlatField)
    referenceFileNames = tomoExperiment.getReferencesFileNames()

    if referenceFileNames == None:
        tomoExperiment.createAverageWfandDf()
        tomoExperiment.findCenterOfRotation()
        print('Cor Found at '+str(tomoExperiment.cor))
        referenceFileNames = tomoExperiment.getReferencesFileNames()

    projectionsFileNames=tomoExperiment.getProjectionsName()
    projectionsFileNames.sort()
    darkFieldFilename=tomoExperiment.getDarkFilename()
    print(darkFieldFilename)
    referenceFileNames.sort()
    print(referenceFileNames)
    ddict={}
    ddict['tomoFileName']=parametersScanFilename
    ddict['projections']=projectionsFileNames
    ddict['references']=referenceFileNames
    ddict['darkField']=darkFieldFilename
    ddict['COR'] = tomoExperiment.cor
    return ddict









if __name__ == "__main__":
    inputFolder='/Volumes/ID17/inhouse4/MITTONE/md1127/AK17-5639/SPECKLES/'
    outputFolder = '/Volumes/ID17/inhouse4/MITTONE/md1127/AK17-5639/SPECKLES/TestOutputSpeckles/'
    mkdir(outputFolder)
    tomoFolders=glob.glob(inputFolder+'/AK175639_speckles_1ms_720p_bin4_bis__003*')
    tomoFolders.sort()
    print(tomoFolders)
    processAllFoldersThreaded(tomoFolders,outputFolder,nbThread=4)

    #result=parseTomoFolderAndCreateRefFiles(tomoFolders[0])
    #print result


