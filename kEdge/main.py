__author__ = 'embrun'
from skimage.feature import register_translation
import EdfFile as edf
import glob
import string
import os
import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk
import pylab
import Image3D





def averageScans(reconstructedFolder,averageFolder):
    folderPath=reconstructedFolder
    listOfFiles=glob.glob(folderPath+'/*.edf')
    listOfFiles.sort()
    numberOfScans=0
    numberOfSlices=0
    minScanNumber = 8192
    for fileName in listOfFiles:
        basename=os.path.basename(fileName)
        splits=basename.split('_')
        scan_number,slice=int(splits[1]),int(splits[3].split('.')[0])
        numberOfScans=max(scan_number,numberOfScans)
        numberOfSlices=max(slice,numberOfSlices)
        minScanNumber=min(scan_number,minScanNumber)

    numberOfScans+=1
    print 'Number Of Scans'+str(numberOfScans)
    listOfScan3D=[]

    for i in range(int(numberOfScans)):
        im=Image3D.Image3D()
        im.scanNumber=i
        im.folderName=folderPath
        im.nbSlices=numberOfSlices
        listOfScan3D.append(im)

    print 'Fin Creation Des Objets'


    for fileName in listOfFiles:
        basename=os.path.basename(fileName)
        splits=basename.split('_')
        scan_number,slice=splits[1],splits[3].split('.')[0]

        listOfScan3D[int(scan_number)].appendFileName(int(slice),fileName)

    print "min Scan Number"+str(minScanNumber)


    for scan in listOfScan3D:
        if not scan.isEmptyScan():
            scan.setDimensions()
            scan.loadSlices()


    scanRef=listOfScan3D[minScanNumber]

    i=0
    for scan in listOfScan3D:
        if(i>minScanNumber):
            scan.stackRegisterWithRef(scanRef)
        i+=1
    i=0
    cpt=0
    for scan in listOfScan3D:
        if(i>minScanNumber and i<numberOfScans-1) and (i !=43)and (i !=44) and (i !=45):
            scanRef.sumImage(scan)
            cpt+=1
        i+=1


    print 'We gonna divide by '+str(cpt)
    scanRef.divide(cpt)
    scanRef.save3DImage(averageFolder+'Reg_Average_'+'_')



if __name__ == "__main__":
    print 'On demarre'
    folderPath='/Volumes/ID17/u836/u836a/helene_juillet2015/AnalysesManu/Below/tomo/EDF_PYHST/'
    resultFolder='/Users/embrun/Essai/Below/Below_'
    averageScans(folderPath,resultFolder)
    folderPath='/Volumes/ID17/u836/u836a/helene_juillet2015/AnalysesManu/Above/tomo/EDF_PYHST/'
    resultFolder='/Users/embrun/Essai/Above/Above_'
    averageScans(folderPath,resultFolder)









