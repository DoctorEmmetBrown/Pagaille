from threading import Thread
from NoiseTracking.OpticalFlow import PavlovSingleMaterial
import os
from InputOutput.pagailleIO import openImage,saveEdf
import numpy as np





class PavlovOpticalFlowSolverThread(Thread):
    def __init__(self,listOfProjections,listOfSampleFileNames,referenceName,darkn,outputFolder):
        Thread.__init__(self)
        self.referenceName = referenceName
        self.listOfProjection=listOfProjections
        self.output=outputFolder
        self.listOfSampleFileNames=listOfSampleFileNames
        self.darkFileName=darkn


    def run(self):
        projectionFiles = self.listOfSampleFileNames
        for numeroProjection in self.listOfProjection:
            numeroProjection = int(numeroProjection)
            print('Processing ' + str(numeroProjection))
            projectionFileName=projectionFiles[numeroProjection]
            result=PavlovSingleMaterial.processOneImage(projectionFileName,self.referenceName,self.darkFileName)
            saveEdf(result, self.output+ '/' + os.path.basename(projectionFileName))

