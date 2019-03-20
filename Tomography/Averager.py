import numpy as np
#from EdfFile import EdfFile as edfFile
from InputOutput.pagailleIO import openImage,saveEdf

import math


class Averager(object):
    def __init__(self, filenames, outputFilename, option=1, parent=None):
        self.fileNames = filenames
        self.option = option
        self.outputFilename = outputFilename
        if (len(self.fileNames) > 0):
            if (option == 1):
                self.processAverage()
            else:
                self.processMedian()
            self.save()

    def processAverage(self):
        print('Average ')
        dat = openImage(self.fileNames[0])
        shape = dat.shape
        self.output = np.zeros(shape)
        for filenam in self.fileNames:
            dat = openImage(filenam)
            self.output += dat
        self.output=np.divide(self.output, len(self.fileNames))

    def processMedian(self):
        print('Median ')

        dat = openImage(self.fileNames[0])
        shape = dat.shape
        self.bigMat = np.zeros((len(self.fileNames), dat.shape[0], dat.shape[1]))
        cpt = 0
        for filenam in self.fileNames:
            dat = openImage(filenam)
            self.bigMat[cpt] += dat
            cpt += 1

        self.output = np.median(self.bigMat, axis=0)

    def save(self):
        print('Save output ' + self.outputFilename)
        saveEdf(self.output,self.outputFilename)


if __name__ == "__main__":
    import sys
    import glob
    folder='/Volumes/BM05/imaging/embrun/160307_gado/Phantom_aboveKedge/'
    filenames=glob.glob(folder+'ref*_3000.edf')
    print(filenames)
    outputfilename='/Volumes/BM05/imaging/embrun/160307_gado/Phantom_aboveKedge/refHST3000.edf'
    av=Averager(filenames,outputfilename,option=0)






