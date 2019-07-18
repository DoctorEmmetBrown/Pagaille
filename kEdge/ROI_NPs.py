import matplotlib.pyplot as plt
from InputOutput.pagailleIO import openImage
import Image3D
import numpy as np
import glob
import os


#5 materials are considered : H2O, I, Ba , Gd, Au,
#10 energies were performed
#                  I-    I+     Ba-   Ba+   Gd-   Gd+  Au-   Au+
mu_water=       [0.341, 0.316,0.294,0.279,0.229,0.224,0.184,0.182]
mu_iodine=      [7.107,33.190,28.200,24.560,12.830,11.550,3.543,3.316]
mu_barium=      [8.217,7.000,5.923,27.320,14.350,12.930,4.000,3.745]
mu_gadolinium=  [12.320,10.490,8.873,7.705,4.015,17.700,5.624,5.272]
mu_gold=        [22.910,19.570,16.590,14.430,7.551,6.809,2.204,8.628]

def registerAndKEdgeSubstractCoupleofImage(belowImage,aboveImage,mum1ebelow,mum2below,mum1above,mum2above):
    print('Registration')

    #belowImage.registerImages(aboveImage)
    denominateur=mum1above*mum2below-mum1ebelow*mum2above
    print('Denominateur:'+str(denominateur))
    concentrationImage=Image3D.Image3D(nbSlices=belowImage.nbSlices,width=belowImage.width,height=belowImage.height)
    concentrationImage.data=np.subtract(mum2below*aboveImage.data,mum2above*belowImage.data)
    concentrationImage.data=1000*concentrationImage.data/denominateur #in mg/mL
    return concentrationImage


def seuilImageFloat(imgFloat, mu):
    print('Seuillage ... ')
    #imgFloat = '/Volumes/ID17/rsrm/md1011/volfloat/Ex_vivo_R0842_16_AuC_PM__AboveAu__001_pag/Ex_vivo_R0842_16_AuC_PM__AboveAu__001_pag_0260.edf'
    img_float = openImage(imgFloat)
    #basename = os.path.basename(imgFloat)
    #fig, ax = plt.subplots(nrows=3, figsize=(120, 20))

    plt.imshow(img_float)
    img_float[img_float<mu]=0
    img_float[img_float>=mu]=1
    plt.imshow(img_float, cmap='gray')

    img_seuil = filtreParTaille(img_float, 4500)
    plt.imshow(img_seuil)
    plt.imshow(img_seuil, cmap='gray')
    #plt.title(basename)

    #plt.show()
    return img_seuil


def loadImages(ddict):
    print('LoadImages')
    numberOfCouples=ddict['nbElements']
    print('numberOfCouples : '+str(numberOfCouples))
    for i in range(0,numberOfCouples):
        element=ddict['listOfElements'][i]

        print('Element:'+str(element))
        coupleOfFolder=ddict['listOfCouplesOfEnergies'][i]

        if len(coupleOfFolder)<2:
            print('PROBLEM!!!!!!!!!')

        else:
            if ('Above' in coupleOfFolder[0]) or (('above' in coupleOfFolder[0])):
                aboveFolder=coupleOfFolder[0]
                basename = os.path.basename(coupleOfFolder[0]).split('Above')[0]
                belowFolder=coupleOfFolder[1]

            else:
                aboveFolder=coupleOfFolder[1]
                belowFolder=coupleOfFolder[0]

            if element=='I' or element=='i' or element=='Iodine' or element=='iodine':
                muM1ebelow=mu_iodine[0]
                muM2ebelow=mu_water[0]
                muM1eabove=mu_iodine[1]
                muM2eabove=mu_water[1]

            if element=='Ba' or element=='ba' or element=='Barium' or element=='barium':
                muM1ebelow=mu_barium[2]
                muM2ebelow=mu_water[2]
                muM1eabove=mu_barium[3]
                muM2eabove=mu_water[3]

            if element=='Gd' or element=='gd' or element=='Gadolinium' or element=='gadolinium':
                muM1ebelow=mu_gadolinium[4]
                muM2ebelow=mu_water[4]
                muM1eabove=mu_gadolinium[5]
                muM2eabove=mu_water[5]

            if element=='Au' or element=='au' or element=='Gold' or element=='gold':
                muM1ebelow=mu_gold[6]
                muM2ebelow=mu_water[6]
                muM1eabove=mu_gold[7]
                muM2eabove=mu_water[7]

        print('AboveFolder = '+str(aboveFolder))
        print('mu Element:'+str(muM1eabove))
        print('mu Water:'+str(muM2eabove))


        aboveImage=Image3D.Image3D(folderName=aboveFolder)
        aboveImagecopy=Image3D.Image3D(folderName=aboveFolder)
        aboveImage.createListOfFiles('tif')
        aboveImage.loadSlices()
        aboveImage.imresize(0.25)


        aboveImage.data=np.asarray(aboveImage.data,np.float32)
        aboveImagecopy.copyImageFromNumpyArray(aboveImage.data)
        print('Above resampling '+str(float(ddict['minValue' + element + 'above']))+' '+str(float(ddict['maxValue' + element + 'above'])))
        aboveImage.resampleIn32Bit(float(ddict['minValue' + element + 'above']),float(ddict['maxValue' + element + 'above']))
        aboveImage.save3DImage('/Users/zaouak/Desktop/Test/AboveFlottant/')
        aboveImage.seuil(0.364)
        aboveImage.filtreParTaille(3000)
        #aboveImage.Dilatation(aboveImage)

        print(' ')
        print('BelowFolder = '+str(belowFolder))
        print('mu Element:'+str(muM1ebelow))
        print('mu Water:'+str(muM2ebelow))
        #outputDirectory=ddict['samplePath']+'/Concentration'+element

        outputDirectory='/Users/zaouak/Desktop/sortie/'
        belowImage=Image3D.Image3D(folderName=belowFolder)
        belowImagecopy=Image3D.Image3D(folderName=belowFolder)

        belowImage.createListOfFiles('tif')
        belowImage.loadSlices()
        belowImage.imresize(0.25)


        belowImage.data=np.asarray(belowImage.data,np.float32)
        belowImagecopy.copyImageFromNumpyArray(belowImage.data)
        print('Below resampling ' + str(float(ddict['minValue' + element + 'below'])) + ' ' + str(float(ddict['maxValue' + element + 'below'])))
        belowImage.resampleIn32Bit(float(ddict['minValue' + element + 'below']),float(ddict['maxValue' + element + 'below']))
        belowImage.save3DImage('/Users/zaouak/Desktop/Test/BelowFlottant/')
        belowImage.seuil(0.364)
        belowImage.filtreParTaille(3000)
        #belowImage.Dilatation(belowImage)


        #belowImagecopy.registerImagesMovingTranslation(belowImage, aboveImage, aboveImagecopy)
        #belowImagecopy.save3DImage('/Users/zaouak/Desktop/Test/BelowRecale/')

        belowImagecopy.registerImagesMovingTranslationAndRotation(belowImage, aboveImage, aboveImagecopy)
        belowImagecopy.save3DImage('/Users/zaouak/Desktop/Test/BelowRecaleRecale/')


        elementConcentrationImage=registerAndKEdgeSubstractCoupleofImage(belowImagecopy,aboveImagecopy,muM1ebelow,muM2ebelow,muM1eabove,muM2eabove)
        elementConcentrationImage.save3DImage('/Users/zaouak/Desktop/Test/Concentration/')
        #elementConcentrationImage.seuilROI(7)

        if(not(os.path.exists(outputDirectory))):
            os.mkdir(outputDirectory)
        elementConcentrationImage.save3DImage(outputDirectory+'/ConcentrationTestRotation_'+ basename + element+'_'+ddict["sampleName"])

        print('---------------')


def browseSampleFolder(folderPath):

    ddict={}
    print('Analysing ... '+str(folderPath))
    sampleName=os.path.basename(folderPath)
    reconstructionsFolder=glob.glob(folderPath + '/Ex_vivo_R0842_16_*')
    reconstructionsFolder.sort()

    constructionFolder = []
    for folder in reconstructionsFolder:
        if ('bin'  not in folder) and ('pag' in folder):
            constructionFolder.append(folder)
    print('construction : ' + str(constructionFolder))

    listOfElements=[]
    ddict['samplePath']=folderPath
    for energyFolder in constructionFolder:
    #    print energyFolder
        basename=os.path.basename(energyFolder)
        print('basename : ' + str(basename))
        #print basename
        if 'Above' in basename:
            elementsOfInterest=(energyFolder.split('Above')[1]).split('_')[0]
            listOfElements.append(elementsOfInterest)
        elif 'above' in basename:
            elementsOfInterest=(energyFolder.split('above')[1]).split('_')[0]
            listOfElements.append(elementsOfInterest)
    #print (listOfElements)

    ddict["sampleName"]=sampleName
    ddict["nbElements"]=len(listOfElements)
    ddict['listOfElements']=listOfElements



    listOfCouplesOfEnergies=[]


    for element in listOfElements:
        #print (element)
        coupleOfEnergyFolder=[]
        for energyFolder in constructionFolder:
            basename=os.path.basename(energyFolder)
            #print('leesgflyegfyueglfl : ' + str(energyFolder))

            fichier = open(energyFolder+ '/Test36/reconstruction_log.info','r')
            Lines = fichier.read().splitlines()

            compteur=0
            for i in range(0,len(Lines)):
                if Lines[i].__contains__("minmax")==True:
                    compteur = i
                    i = len(Lines)

            values = Lines[compteur]
            #print(values)


            values=values.split('[')
            values=values[1]
            values=values.split(',')
            values1=values[1]
            values1=values1.split(']')
            values1=values1[0]

            values2=[]
            values2.append(values[0])
            values2.append(values1)
            #print(values2)

            minValue = values2[0]
            maxValue = values2[1]

            print('----------------')
            #print(sampleName)
            print('minValue : '+str(minValue))
            print('maxValue : '+str(maxValue))
            print('----------------')

            if 'Above'+element in basename:
                coupleOfEnergyFolder.append(energyFolder + '/Test36')
                ddict['minValue' + element + 'above'] = float(minValue)
                ddict['maxValue' + element + 'above'] = float(maxValue)
            if 'above'+element in basename:
                coupleOfEnergyFolder.append(energyFolder + '/Test36')
                ddict['minValue' + element + 'above'] = float(minValue)
                ddict['maxValue' + element + 'above'] = float(maxValue)
            if 'below'+element in basename:
                coupleOfEnergyFolder.append(energyFolder + '/Test36')
                ddict['minValue' + element + 'below'] = float(minValue)
                ddict['maxValue' + element + 'below'] = float(maxValue)
            if 'Below'+element in basename:
                coupleOfEnergyFolder.append(energyFolder + '/Test36')
                ddict['minValue' + element + 'below'] = float(minValue)
                ddict['maxValue' + element + 'below'] = float(maxValue)
            #print(coupleOfEnergyFolder)

            fichier.close()
        #print coupleOfEnergyFolder
        listOfCouplesOfEnergies.append(coupleOfEnergyFolder)

    #print listOfCouplesOfEnergies
    ddict['listOfElements'] = listOfElements
    ddict['listOfCouplesOfEnergies']=listOfCouplesOfEnergies

    return ddict


if __name__ == "__main__":
    print('Seuillage + Recalage')
    print('------------')

    secondfolder = '/Users/zaouak/Desktop/Test'
                   #'/Volumes/ID17/rsrm/md1011/voltif/'
    secondsamplesList = glob.glob(secondfolder)
    secondsamplesList.sort()

    for sample in secondsamplesList:
        dictForElement = browseSampleFolder(sample)
        print(dictForElement)
        loadImages(dictForElement)
        print('-------------------')
        print('-------------------')
        print('\n')