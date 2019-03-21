# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:26:20 2016

@author: broche

"""

import os
import numpy as np
import PyMcaQt as qt
import SimpleITK as sitk
import pyqtgraph as pg
import time
import copy
from skimage import transform
from scipy.interpolate import splprep, splev
from Interactor3D import Interactor3D
from ImportThread import ImportThread, ImportNo, ImportDicom
from ExportThread import ExportThread
from LabelEditAndButton import LabelEditAndButton
from OutputDialog import OutputDialog
from SliderAndLabel import SliderAndLabel
from ImageReader import DicomReader, MatReader, STLReader
from RegistrationGUI import RegisteringOption
from Reader_GUI import DicomReaderGUI, MatReaderGUI
from VolumeRenderingGUI import VolumeRenderingGUI
import ImageProcessing as IP
import RegisteringIP as IPR
import registerThread as rT
import glob
import nibabel as nib
from scipy.cluster.vq import kmeans2, whiten
import scipy.stats as stat
import ImageWritter


class Frame(qt.QWidget):

    def __init__(self, parent=None):

        qt.QWidget.__init__(self, parent)

        """ Attributs """
        self.parent = parent
        self.main_path = '/data/id17/broncho/'
        self.loaded_path = '/'

        self.Data_list = []
        self.Name_list = []

        self.ItemsLists = []
        self.ItemsInit = {'Seeds': {'Direction0': [], 'Direction1': [], 'Direction2': []}, \
                          'Zones': {'Direction0': [], 'Direction1': [], 'Direction2': []}, \
                          'Circles': {'Direction0': [], 'Direction1': [], 'Direction2': []}, \
                          'Arrows': {'Direction0': [], 'Direction1': [], 'Direction2': []}, \
                          'Poly': {'Direction0': [[]], 'Direction1': [[]], 'Direction2': [[]]}}

        self.Overlays = []
        self.OverlayPar = {'Flag': False, 'Range': [-1, -1], 'Alpha': 0.5, 'Image': [], 'ColorMap': []}

        self.listDicomHeader = {}

        """ Widgets Initialisation """

        self.tabWidget = qt.QTabWidget()
        self.tabPlanesView = qt.QWidget()
        self.tab3DViewer = qt.QWidget()

        self.plt = pg.GraphicsWindow()
        self.plt.hide()

        self.mainLayout = qt.QGridLayout()

        self.layoutViewPlanes = qt.QGridLayout()
        self.buttonLayout = qt.QVBoxLayout()
        self.image3DWidget = Interactor3D(self)
        self.imageSelection = qt.QListWidget(self)
        self.imageSelection.setContextMenuPolicy(qt.Qt.CustomContextMenu)
        # self.connect(self.imageSelection, qt.SIGNAL("customContextMenuRequested(QPoint)"),  self.listItemRightClicked)
        self.imageSelection.customContextMenuRequested.connect(self.listItemRightClicked)
        self.progressBar = qt.QProgressBar(self)
        self.mouseDisplay = LabelEditAndButton(True, "", False, "", False, "")

        self.layout3d = qt.QGridLayout()
        self.VRGUI = VolumeRenderingGUI()

        """ Widgets Design"""

        self.imageSelection.setMaximumWidth(250)

        """Signals"""

        self.imageSelection.currentRowChanged.connect(self._dataToShowChanged)
        self.image3DWidget.axialWidget.MovedOnVizualizer.connect(self._moved)
        self.image3DWidget.coronalWidget.MovedOnVizualizer.connect(self._moved)
        self.image3DWidget.sagittalWidget.MovedOnVizualizer.connect(self._moved)

        self.image3DWidget.axialWidget.clickedOnVizualizer.connect(self._cliked)
        self.image3DWidget.coronalWidget.clickedOnVizualizer.connect(self._cliked)
        self.image3DWidget.sagittalWidget.clickedOnVizualizer.connect(self._cliked)

        self.image3DWidget.axialWidget.releasedOnVizualizer.connect(self._released)
        self.image3DWidget.coronalWidget.releasedOnVizualizer.connect(self._released)
        self.image3DWidget.sagittalWidget.releasedOnVizualizer.connect(self._released)

        self.image3DWidget.toolBar.radius.lineEdit.textChanged.connect(self._changeRadius)

        self.image3DWidget.toolBar.pointRemoveAction.triggered.connect(self._RemoveItems)

        self.tabWidget.currentChanged.connect(self._tabChanged)

        """Placing widgets"""

        self.buttonLayout.addWidget(self.imageSelection)
        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.layoutViewPlanes.addWidget(self.image3DWidget, 0, 0)
        self.layoutViewPlanes.addWidget(self.plt, 0, 1)
        self.layoutViewPlanes.addLayout(self.buttonLayout, 0, 2)
        self.layoutViewPlanes.addWidget(self.progressBar, 1, 0)
        self.layoutViewPlanes.addWidget(self.mouseDisplay, 1, 2)
        self.tabPlanesView.setLayout(self.layoutViewPlanes)

        self.layout3d.addWidget(self.VRGUI)
        self.tab3DViewer.setLayout(self.layout3d)

        self.tabWidget.addTab(self.tabPlanesView, 'Slices Viewer')
        self.tabWidget.addTab(self.tab3DViewer, 'Volumetric Viewer')

        self.mainLayout.addWidget(self.tabWidget)

        self.setLayout(self.mainLayout)

    """ Menu Methods """


    """Internal Methods"""

    def listItemRightClicked(self, QPos):
        self.listMenu = qt.QMenu()
        remove_item = self.listMenu.addAction("Delete")
        remove_item.triggered.connect(self.removeImage)
        parentPosition = self.imageSelection.mapToGlobal(qt.QPoint(0, 0))
        self.listMenu.move(parentPosition + QPos)
        self.listMenu.show()

    def removeImage(self):

        del self.Data_list[self.imageSelection.currentRow()]
        del self.Name_list[self.imageSelection.currentRow()]
        del self.ItemsLists[self.imageSelection.currentRow()]
        del self.Overlays[self.imageSelection.currentRow()]

        self.imageSelection.takeItem(self.imageSelection.currentRow())
        self.imageSelection.setCurrentRow(self.imageSelection.count() - 1)
        self.setLayout(self.mainLayout)

    """I/O Methods"""

    def _load(self, inputFolder=-1):

        if inputFolder == False:
            self.inputFiles = qt.QFileDialog.getOpenFileNames(self, "Select one or more files to open", self.main_path)
            self.inputFiles = self.inputFiles[0]
            if not self.inputFiles:
                return 0
        else:

            inputFolder += "/Paganin/EDF/"
            self.inputFiles = [f for f in os.listdir(inputFolder) if os.path.isfile(os.path.join(inputFolder, f))]
            self.inputFiles.sort()

            for i in range(0, len(self.inputFiles)):
                self.inputFiles[i] = inputFolder + '/' + self.inputFiles[i]

        if len(self.inputFiles) != 0:
            self.Name_list.append(os.path.basename(self.inputFiles[0]))

        path, filename = os.path.split(self.inputFiles[0][0])

        self.progressBar.setRange(0, len(self.inputFiles) - 1)
        self.importThread = ImportThread(self.inputFiles, self)
        self.importThread.Progress.connect(self.setProgressBar)
        self.importThread.finished.connect(self.setData)

        self.importThread.start()

    def _loadNoThread(self, inputFolder=-1):
        if inputFolder == -1:
            self.inputFiles = qt.QFileDialog.getOpenFileNames(self, "Select one or more files to open", self.main_path)
        else:

            # inputFolder += '/Paganin/EDF/'
            self.inputFiles = [f for f in os.listdir(inputFolder) if os.path.isfile(os.path.join(inputFolder, f))]
            self.inputFiles.sort()

            for i in range(0, len(self.inputFiles)):
                self.inputFiles[i] = inputFolder + '/' + self.inputFiles[i]

        self.Name_list.append(self.inputFiles[0].split('/')[-1])
        self.importedData = ImportNo(self.inputFiles)
        self.setDataNo()

    def _loadFolders(self):
        self.inputDirectory = str(qt.QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.inputDirectories = [d for d in os.listdir(self.inputDirectory) if
                                 os.path.isdir(os.path.join(self.inputDirectory, d))]

        self.inputDirectories.sort()
        self.foldersToOpen = DicomReaderGUI(self.inputDirectories)

        self.foldersToOpen.startImporting.clicked.connect(self._startLoadingFolders)
        self.foldersToOpen.show()

    def _startLoadingFolders(self):

        self.foldersToOpen.hide()
        for w in self.foldersToOpen.widgetList:
            if w.checkState() == 2:
                try:
                    self._loadNoThread(self.inputDirectory + '/' + w.text() + '/')
                except:
                    print('Not Imported')

    def _loadSTL(self):
        self.inputFile = qt.QFileDialog.getOpenFileNames(self, "Select STL File to Open", self.main_path)

        reader = STLReader(self.inputFile[0])
        stlData = reader.data()
        self.Name_list.append('Mesh')

        self.ItemsLists.append(copy.deepcopy(self.ItemsInit))
        self.Overlays.append(copy.deepcopy(self.OverlayPar))

        self.Data_list.append(stlData)

        self.image3DWidget._setDataVolume(self.Data_list[self.imageSelection.currentRow()])
        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])
        self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])
        self.progressBar.reset()

        item = qt.QListWidgetItem(self.Name_list[-1])
        item.setTextAlignment(qt.Qt.AlignHCenter)
        item.setFlags(item.flags() | qt.Qt.ItemIsEditable)
        self.imageSelection.itemChanged.connect(self.nameImageChange)
        self.imageSelection.addItem(item)

        self.imageSelection.setMaximumWidth(250)
        self.imageSelection.setCurrentRow(self.imageSelection.count() - 1)
        self.setLayout(self.mainLayout)

    def _loadDicom(self):

        self.inputFiles = qt.QFileDialog.getOpenFileNames(self, "Select Dicom/Mat Files to Open", self.main_path)
        print('*-*-')
        print(self.inputFiles[0])


        self.loaded_path = os.path.dirname(self.inputFiles[0][0])

        if self.inputFiles[0][0].endswith('.mat'):
            self.matClass = MatReader(self.inputFiles[0][0])
            self.matGUI = MatReaderGUI(self.matClass.getListScan())
            self.buttonLayout.addWidget(self.matGUI)
            self.matGUI.startImporting.clicked.connect(self._startLoadingMat)
        elif self.inputFiles[0][0].endswith('.nrrd'):

            data, xxxx = nrrd.read(self.inputFiles[0][0])
            data = np.array(data)

            self.addImage(self.inputFiles[0][0], np.swapaxes(data, 0, 1))

        elif self.inputFiles[0][0].endswith('.nii'):

            filenii = nib.load(self.inputFiles[0][0])
            ImageOut = np.asarray(filenii.get_data())
            shapeIm = ImageOut.shape

            if len(shapeIm) > 3:
                nImage = shapeIm[3]
                for ni in range(0, nImage):
                    self.addImage(str(ni) + '_' + self.inputFiles[0][0].split('/')[-1], ImageOut[:, :, :, :, ni])

            elif len(shapeIm) <= 3:
                self.addImage(self.inputFiles[0][0].split('/')[-1], ImageOut)




        else:
            quit_msg = "Are the images ordered?"
            reply = qt.QMessageBox.question(self, 'Message', quit_msg, qt.QMessageBox.Yes, qt.QMessageBox.No)

            if reply == qt.QMessageBox.Yes:

                self.dicomClass = DicomReader([self.inputFiles[0][0]])
                self.dicomClass.getListScan()
                a = self.dicomClass.a
                b = self.dicomClass.b
                self.patientName = self.dicomClass.patientName
                self.parent.setWindowTitle(self.patientName)

                self.importedData = ImportDicom(self.inputFiles).inputData
                px_z = self.dicomClass.pixel_size[0]
                px_x = self.dicomClass.pixel_size[1]
                px_y = self.dicomClass.pixel_size[2]

                factor = [1.0 / px_x, 1.0 / px_y, 1.0 / px_z]

                imageOut = a * IP.resizeImage(np.copy(self.importedData), "Linear", factor) + b
                self.addImage(self.dicomClass.SerieDescription, imageOut, "")
                self.image3DWidget._setDataVolume(imageOut, -150, 300)

            else:
                self.dicomClass = DicomReader(self.inputFiles)

                if not len(self.dicomClass.getListScan()) == 0:
                    self.dicomGUI = DicomReaderGUI(self.dicomClass.getListScan())
                    self.patientName = self.dicomClass.patientName
                    self.parent.setWindowTitle(self.patientName)

                    self.dicomGUI.startImporting.clicked.connect(self._startLoadingDicom)
                    self.buttonLayout.addWidget(self.dicomGUI)


                else:
                    self.dicomClass.getListScan()
                    self.patientName = self.dicomClass.patientName
                    self.parent.setWindowTitle(self.patientName)
                    dataToImport = self.dicomClass.data
                    self.addImage(self.inputFiles[0], dataToImport)

    def _startLoadingDicom(self):
        self.dicomGUI.hide()
        self.listToImport = []

        for w in self.dicomGUI.widgetList:
            if w.checkState() == 2:
                self.listToImport.append(w.text())

        dataToImport = self.dicomClass.importScan(self.listToImport)
        self.listDicomHeader.update(self.dicomClass.info)

        for i, scan in enumerate(dataToImport):
            px_z = self.dicomClass.pixel_size[0]
            px_x = self.dicomClass.pixel_size[1]
            px_y = self.dicomClass.pixel_size[2]

            factor = [1.0 / px_x, 1.0 / px_y, 1.0 / px_z]

            imageOut = IP.resizeImage(np.copy(dataToImport[scan]), "Linear", factor)
            self.addImage(scan, imageOut, str(self.listDicomHeader[scan]))

            self.image3DWidget._setDataVolume(imageOut, -150, 300)

        del dataToImport

    def _startLoadingMat(self):

        self.matGUI.hide()
        self.listToImport = []

        for w in self.matGUI.widgetList:
            if w.checkState() == 2:
                self.listToImport.append(w.text())

        for scan in self.listToImport:
            data = self.matClass.dicmat[scan]

            if data.ndim < 3:

                image = np.zeros((1, data.shape[0], data.shape[1]))
                image[0, :, :] = data
                self.addImage(scan, image)
            else:
                self.addImage(scan, self.matClass.dicmat[scan])

    def _save(self):
        dialog = OutputDialog(self)
        ret = dialog.exec_()

        if (ret == qt.QDialog.Accepted):
            extensionNumber = dialog.outputImageFormat.currentIndex()
            if (extensionNumber == 0):
                extension = '.edf'
            elif (extensionNumber == 1):
                extension = '.tif'
            elif (extensionNumber == 2):
                extension = '.nii'
            elif (extensionNumber == 3):
                extension = '.npy'
            elif (extensionNumber == 4):
                extension = '.mat'
            elif (extensionNumber == 5):
                extension = '.png'
            elif (extensionNumber == 6):
                extension = '.dcm'
            elif (extensionNumber == 7):
                extension = '.nrrd'

            self.resultFileName = qt.QFileDialog.getSaveFileName(self, "Save Image Sequence ", self.loaded_path,
                                                                 'Images (*' + extension + ')')[0]

            dataToStore = self.Data_list[self.imageSelection.currentRow()]

            self.progressBar.setRange(0, dataToStore.shape[0])

            Z = self.Data_list[self.imageSelection.currentRow()].shape[0]
            X = self.Data_list[self.imageSelection.currentRow()].shape[1]
            Y = self.Data_list[self.imageSelection.currentRow()].shape[2]

            if extensionNumber == 2:

                ImageWritter.write_nifti(np.swapaxes(dataToStore, 0, 2), self.resultFileName + '.nii')
            elif extensionNumber == 7:

                nrrd.write(self.resultFileName + '.nrrd', np.swapaxes(dataToStore, 0, 2))
            else:
                self.exportThread = ExportThread(self.resultFileName, dataToStore, extension, self)
                self.exportThread.Progress.connect(self.setProgressBar)
                self.exportThread.start()

    """ Signals Methods"""

    def _tabChanged(self, tabIndex):
        if tabIndex == 1:
            self.VRGUI.ImagesList = self.Name_list
            self.VRGUI.setImages()
            self.VRGUI.DataList = self.Data_list
            self.VRGUI.ItemLists = self.ItemsLists

    def _changeRadius(self, radius):
        newRadius = float(radius)

        if self.image3DWidget.toolBar.drawingAction.isChecked() == True:

            if self.circleToMove[0] == 0:
                self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction0'][-1][3] = newRadius
            if self.circleToMove[0] == 1:
                self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction1'][-1][3] = newRadius
            if self.circleToMove[0] == 2:
                self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction2'][-1][3] = newRadius

            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

    def setData(self):

        self.ItemsLists.append(copy.deepcopy(self.ItemsInit))
        self.Overlays.append(copy.deepcopy(self.OverlayPar))

        self.Data_list.append(self.importThread.inputData)

        self.image3DWidget._setDataVolume(self.Data_list[self.imageSelection.currentRow()])
        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])
        self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])
        self.progressBar.reset()

        item = qt.QListWidgetItem(self.Name_list[-1])
        item.setTextAlignment(qt.Qt.AlignHCenter)
        item.setFlags(item.flags() | qt.Qt.ItemIsEditable)
        self.imageSelection.itemChanged.connect(self.nameImageChange)
        self.imageSelection.addItem(item)

        self.imageSelection.setMaximumWidth(250)
        self.imageSelection.setCurrentRow(self.imageSelection.count() - 1)
        self.setLayout(self.mainLayout)

    def setDataNo(self):

        self.ItemsLists.append(copy.deepcopy(self.ItemsInit))
        self.Overlays.append(copy.deepcopy(self.OverlayPar))

        self.Data_list.append(self.importedData.inputData)

        self.image3DWidget._setDataVolume(self.Data_list[self.imageSelection.currentRow()])
        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])
        self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])
        self.progressBar.reset()

        item = qt.QListWidgetItem(self.Name_list[-1])
        item.setTextAlignment(qt.Qt.AlignHCenter)
        item.setFlags(item.flags() | qt.Qt.ItemIsEditable)
        self.imageSelection.itemChanged.connect(self.nameImageChange)
        self.imageSelection.addItem(item)

        self.imageSelection.setMaximumWidth(250)
        self.imageSelection.setCurrentRow(self.imageSelection.count() - 1)
        self.setLayout(self.mainLayout)

    def nameImageChange(self):
        self.Name_list[self.imageSelection.currentRow()] = self.imageSelection.currentItem().text()

    def setProgressBar(self, i):
        self.progressBar.setValue(i)

    def _dataToShowChanged(self):

        self.image3DWidget._setDataVolume(self.Data_list[self.imageSelection.currentRow()])
        if len(self.ItemsLists) > self.imageSelection.currentRow():
            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

        if len(self.Overlays) > self.imageSelection.currentRow():
            self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])

    def _moved(self, ddict):

        x = ddict['x']
        y = ddict['y']
        z = ddict['z']

        try:
            stringMouse = '(' + str(x) + ' , ' + str(y) + ' , ' + str(z) + ')     '
            textValue = '%.3f' % (self.Data_list[self.imageSelection.currentRow()][z, y, x])
            stringMouse += textValue

        except:
            stringMouse = '(' + str(int(ddict['x'])) + ' , ' + str(int(ddict['y'])) + ' , ' + str(z) + ')  ------'

        self.mouseDisplay.changeLabel(stringMouse)

    def _cliked(self, ddict):
        x = ddict['x']
        y = ddict['y']
        z = ddict['z']
        PlaneSection = ddict['PlaneSection']
        if self.image3DWidget.toolBar.pointerAction.isChecked() == True:
            seed = [x, y, z]

            self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction0'].append(seed)
            self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction1'].append(seed)
            self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction2'].append(seed)
            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

        if self.image3DWidget.toolBar.zone1Action.isChecked() == True:
            self.clic_zone = [x, y, z]

        if self.image3DWidget.toolBar.polygonAction.isChecked() == True:

            seed = [x, y, z]
            if (ddict['event'] != 'RMousePressed'):
                if PlaneSection == 0:

                    if len(self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction0'][-1]) != 0:
                        if self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction0'][-1][-1][2] != z:
                            self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction0'].append([])

                        self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction0'][-1].append(
                            [seed[0] + 0.5, seed[1] + 0.5, seed[2]])
                    else:
                        self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction0'][-1].append(
                            [seed[0] + 0.5, seed[1] + 0.5, seed[2]])

                if PlaneSection == 1:
                    if len(self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction1'][-1]) != 0:
                        if self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction1'][-1][-1][1] != y:
                            self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction1'].append([])
                        self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction1'][-1].append(
                            [seed[0] + 0.5, seed[1], seed[2] + 0.5])


                    else:
                        self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction1'][-1].append(
                            [seed[0] + 0.5, seed[1], seed[2] + 0.5])

                if PlaneSection == 2:
                    if len(self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction2'][-1]) != 0:
                        if self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction2'][-1][-1][0] != x:
                            self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction2'].append([])

                        self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction2'][-1].append(
                            [seed[0], seed[1] + 0.5, seed[2] + 0.5])
                    else:
                        self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction2'][-1].append(
                            [seed[0], seed[1] + 0.5, seed[2] + 0.5])
            else:
                if PlaneSection == 0:
                    self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction0'].append([])
                if PlaneSection == 1:
                    self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction1'].append([])
                if PlaneSection == 2:
                    self.ItemsLists[self.imageSelection.currentRow()]['Poly']['Direction2'].append([])

            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

        if self.image3DWidget.toolBar.drawingAction.isChecked() == True:
            seed = [x, y, z, float(self.image3DWidget.toolBar.radius.lineEdit.text())]
            if PlaneSection == 0:
                flag_in_circle = False
                for i, circle in enumerate(self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction0']):
                    if (((circle[0] - x) ** 2 + (circle[1] - y) ** 2) ** 0.5) < (circle[3]):
                        flag_in_circle = True
                        self.circleToMove = [i, x, y]

                if not flag_in_circle:
                    self.circleToMove = [0]
                    self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction0'].append(seed)
            if PlaneSection == 1:
                flag_in_circle = False
                for i, circle in enumerate(self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction1']):
                    if (((circle[0] - x) ** 2 + (circle[2] - z) ** 2) ** 0.5) < (circle[3]):
                        flag_in_circle = True
                        self.circleToMove = [i, x, z]
                if not flag_in_circle:
                    self.circleToMove = [1]
                    self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction1'].append(seed)
            if PlaneSection == 2:
                flag_in_circle = False
                for i, circle in enumerate(self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction2']):
                    if (((circle[1] - y) ** 2 + (circle[2] - z) ** 2) ** 0.5) < (circle[3]):
                        flag_in_circle = True
                        self.circleToMove = [i, y, z]
                if not flag_in_circle:
                    self.circleToMove = [2]
                    self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction2'].append(seed)

            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

    def _released(self, ddict):
        x = ddict['x']
        y = ddict['y']
        z = ddict['z']
        PlaneSection = ddict['PlaneSection']

        if self.image3DWidget.toolBar.zone1Action.isChecked() == True:

            if PlaneSection == 0:
                self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction0'].append(
                    [self.clic_zone[0], self.clic_zone[1], self.clic_zone[2], x, y, z])
            if PlaneSection == 1:
                self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction1'].append(
                    [self.clic_zone[0], self.clic_zone[1], self.clic_zone[2], x, y, z])
            if PlaneSection == 2:
                self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction2'].append(
                    [self.clic_zone[0], self.clic_zone[1], self.clic_zone[2], x, y, z])

            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

        if self.image3DWidget.toolBar.drawingAction.isChecked() == True:
            if PlaneSection == 0:
                if len(self.circleToMove) > 1:
                    circleToMove = self.circleToMove[0]

                    for i, circle in enumerate(
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction0']):
                        if i == circleToMove:
                            dx = x - self.circleToMove[1]
                            dy = y - self.circleToMove[2]
                            self.circleToMove = [0]
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction0'][i][0] += dx
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction0'][i][1] += dy
                            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

            if PlaneSection == 1:
                if len(self.circleToMove) > 1:
                    circleToMove = self.circleToMove[0]

                    for i, circle in enumerate(
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction1']):
                        if i == circleToMove:
                            dx = x - self.circleToMove[1]
                            dz = z - self.circleToMove[2]
                            self.circleToMove = [1]
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction1'][i][0] += dx
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction1'][i][2] += dz
                            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

            if PlaneSection == 2:
                if len(self.circleToMove) > 1:
                    circleToMove = self.circleToMove[0]

                    for i, circle in enumerate(
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction2']):
                        if i == circleToMove:
                            dy = y - self.circleToMove[1]
                            dz = z - self.circleToMove[2]
                            self.circleToMove = [2]
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction2'][i][1] += dy
                            self.ItemsLists[self.imageSelection.currentRow()]['Circles']['Direction2'][i][2] += dz
                            self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

    def _RemoveItems(self):

        self.ItemsLists[self.imageSelection.currentRow()] = {}
        self.ItemsLists[self.imageSelection.currentRow()] = copy.deepcopy(self.ItemsInit)
        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])
        self.image3DWidget.axialWidget._changeSlice()
        self.image3DWidget.coronalWidget._changeSlice()
        self.image3DWidget.sagittalWidget._changeSlice()

    def _constante1(self):
        if self.Image1.currentText() == 'Constant':
            self.constante1.show()
        else:
            self.constante1.hide()

    def _constante2(self):
        if self.Image2.currentText() == 'Constant':
            self.constante2.show()
        else:
            self.constante2.hide()

    def _radius(self):
        if self.InitialLS.currentText() == 'Starting Seed Points':
            self.radius.show()
        else:
            self.radius.hide()

    def Ero(self):
        if self.checkErosion.checkState() == 2:
            self.checkDilatation.setCheckState(0)
            self.checkOpening.setCheckState(0)
            self.checkClosing.setCheckState(0)
            self.checkFilling.setCheckState(0)
            self.checkDistance.setCheckState(0)
            self.checkSkeletton.setCheckState(0)

    def Dil(self):
        if self.checkDilatation.checkState() == 2:
            self.checkErosion.setCheckState(0)
            self.checkOpening.setCheckState(0)
            self.checkClosing.setCheckState(0)
            self.checkFilling.setCheckState(0)
            self.checkDistance.setCheckState(0)
            self.checkSkeletton.setCheckState(0)

    def Open(self):
        if self.checkOpening.checkState() == 2:
            self.checkErosion.setCheckState(0)
            self.checkDilatation.setCheckState(0)
            self.checkClosing.setCheckState(0)
            self.checkFilling.setCheckState(0)
            self.checkDistance.setCheckState(0)
            self.checkSkeletton.setCheckState(0)

    def Closing(self):
        if self.checkClosing.checkState() == 2:
            self.checkErosion.setCheckState(0)
            self.checkDilatation.setCheckState(0)
            self.checkOpening.setCheckState(0)
            self.checkFilling.setCheckState(0)
            self.checkDistance.setCheckState(0)
            self.checkSkeletton.setCheckState(0)

    def Filling(self):
        if self.checkFilling.checkState() == 2:
            self.checkErosion.setCheckState(0)
            self.checkDilatation.setCheckState(0)
            self.checkOpening.setCheckState(0)
            self.checkClosing.setCheckState(0)
            self.checkDistance.setCheckState(0)
            self.checkSkeletton.setCheckState(0)

    def Distance(self):
        if self.checkDistance.checkState() == 2:
            self.checkErosion.setCheckState(0)
            self.checkDilatation.setCheckState(0)
            self.checkOpening.setCheckState(0)
            self.checkClosing.setCheckState(0)
            self.checkSkeletton.setCheckState(0)

    def Skeletton(self):
        if self.checkSkeletton.checkState() == 2:
            self.checkErosion.setCheckState(0)
            self.checkDilatation.setCheckState(0)
            self.checkOpening.setCheckState(0)
            self.checkClosing.setCheckState(0)
            self.checkDistance.setCheckState(0)

    """ Image Processing Function  """

    """ Selection"""

    def _smoothContours(self):

        for i, curve in enumerate(self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction0"]):
            if len(curve) != 0:
                new_arr = []
                plane = curve[0][2]
                for coord in curve:
                    if len(coord) != 0:
                        new_arr.append(coord[0:2])

                tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
                u_new = np.linspace(u.min(), u.max(), len(coord) * 1000)
                x, y = splev(u_new, tck, der=0)

                self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction0"][i] = []
                for j, xj in enumerate(x):
                    self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction0"][i].append([xj, y[j], plane])

        for i, curve in enumerate(self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction1"]):
            if len(curve) != 0:
                new_arr = []
                plane = curve[0][1]
                for coord in curve:
                    if len(coord) != 0:
                        new_arr.append([coord[0], coord[2]])

                tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
                u_new = np.linspace(u.min(), u.max(), len(coord) * 100)
                x, y = splev(u_new, tck, der=0)

                self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction1"][i] = []
                for j, xj in enumerate(x):
                    self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction1"][i].append([xj, plane, y[j]])

        for i, curve in enumerate(self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction2"]):
            if len(curve) != 0:
                new_arr = []
                plane = curve[0][0]
                for coord in curve:
                    if len(coord) != 0:
                        new_arr.append([coord[1], coord[2]])

                tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
                u_new = np.linspace(u.min(), u.max(), len(coord) * 100)
                x, y = splev(u_new, tck, der=0)

                self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction2"][i] = []
                for j, xj in enumerate(x):
                    self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction2"][i].append([plane, xj, y[j]])

        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])

    def _interpolateContourGUI(self):

        self.hide_all_button()
        self.smoothButton = qt.QPushButton("Interpolating")

        self.smoothButton.clicked.connect(self._interpolateContour)

        self.smoothButton.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.smoothButton)

    def _interpolateContour(self):

        contours = self.ItemsLists[self.imageSelection.currentRow()]["Poly"]

        nz = self.Data_list[self.imageSelection.currentRow()].shape[0]
        nx = self.Data_list[self.imageSelection.currentRow()].shape[1]
        ny = self.Data_list[self.imageSelection.currentRow()].shape[2]

        new_Image = IP.InterpolateDataPoints(contours, [nz, nx, ny])

        self.addImage("Interpolated_Contour_", new_Image)

    def _interpolateMaskGUI(self):

        self.hide_all_button()
        self.smoothButton = qt.QPushButton("Smooth")
        self.nbIter = LabelEditAndButton(True, "Smoothing Number of Iterations : ", True, str(20), False)
        self.InterF = LabelEditAndButton(True, "OutputImage Interpolation Factor: ", True, str(4), False)

        self.smoothButton.clicked.connect(self._interpolateMask)

        self.smoothButton.setMaximumWidth(250)
        self.nbIter.setMaximumWidth(250)
        self.InterF.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.nbIter)
        self.buttonLayout.addWidget(self.InterF)
        self.buttonLayout.addWidget(self.smoothButton)

        self.setLayout(self.mainLayout)

    def _interpolateMask(self):

        Image = self.Data_list[self.imageSelection.currentRow()]
        NbIter = int(str(self.nbIter.lineEdit.text()))
        Factor = int(str(self.InterF.lineEdit.text()))
        ImageReturn = IP.smooth_3D(np.copy(Image), NbIter, Factor)
        self.addImage("Smooth_Mask", ImageReturn)

    """ Preprocess CT """

    def correct_ring1_GUI(self):
        self.hide_all_button()
        self.alpha = LabelEditAndButton(True, "Alpha : ", True, str(3.0), False)
        self.block = LabelEditAndButton(True, "Block : ", True, str(2.0), False)
        self.filterButton = qt.QPushButton("Filter")
        self.filterButton.clicked.connect(self.correct_ring1)

        self.alpha.setMaximumWidth(250)
        self.block.setMaximumWidth(250)
        self.filterButton.setMaximumWidth(250)
        self.buttonLayout.addWidget(self.alpha)
        self.buttonLayout.addWidget(self.block)
        self.buttonLayout.addWidget(self.filterButton)
        self.setLayout(self.mainLayout)

    """ Filters """

    def add_alpha_map_GUI(self):

        self.hide_all_button()
        self.addAlphaMap = qt.QPushButton("Add Overlay Map")

        self.alphaValue = LabelEditAndButton(True, "Alpha Value :", True, str(0.5), False)

        self.minValue = LabelEditAndButton(True, "Min Value :", True, str(-1.0), False)
        self.maxValue = LabelEditAndButton(True, "Max Value :", True, str(-1.0), False)

        self.Image2Map = qt.QComboBox()
        self.Image2Map.addItems(self.Name_list)

        self.ColorMapB = qt.QComboBox()
        self.ColorMapB.addItem('GrayLevel')
        self.ColorMapB.addItem('Jet')

        self.addAlphaMap.clicked.connect(self.add_alpha_map)

        self.minValue.setMaximumWidth(250)
        self.maxValue.setMaximumWidth(250)
        self.addAlphaMap.setMaximumWidth(250)
        self.Image2Map.setMaximumWidth(250)
        self.ColorMapB.setMaximumWidth(250)
        self.alphaValue.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.Image2Map)
        self.buttonLayout.addWidget(self.ColorMapB)
        self.buttonLayout.addWidget(self.alphaValue)
        self.buttonLayout.addWidget(self.minValue)
        self.buttonLayout.addWidget(self.maxValue)
        self.buttonLayout.addWidget(self.addAlphaMap)
        self.setLayout(self.mainLayout)

    def add_alpha_map(self):

        ImageToMap = self.Data_list[self.Image2Map.currentIndex()]
        AlphaValue = float(str(self.alphaValue.lineEdit.text()))
        colorMap = self.image3DWidget.toolBar.colormapList[self.ColorMapB.currentIndex()]

        minV = float(str(self.minValue.lineEdit.text()))
        maxV = float(str(self.maxValue.lineEdit.text()))

        self.Overlays[self.imageSelection.currentRow()]["Flag"] = True
        self.Overlays[self.imageSelection.currentRow()]["Range"][0] = minV
        self.Overlays[self.imageSelection.currentRow()]["Range"][1] = maxV
        self.Overlays[self.imageSelection.currentRow()]["Alpha"] = AlphaValue
        self.Overlays[self.imageSelection.currentRow()]["Image"] = ImageToMap
        self.Overlays[self.imageSelection.currentRow()]["ColorMap"] = colorMap

        self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])

    def remove_alpha_map(self):
        self.Overlays[self.imageSelection.currentRow()] = copy.deepcopy(self.OverlayPar)
        self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])

    def _histoGUI(self):
        self.hide_all_button()
        self.histoButton = qt.QPushButton("Histo")
        self.minBin = LabelEditAndButton(True, "Min Bin Value :", True, str(0), False)
        self.maxBin = LabelEditAndButton(True, "Max Bin Value :", True, str(40.0), False)
        self.nbBin = LabelEditAndButton(True, "Number of Bin :", True, str(10000), False)
        self.percentile = LabelEditAndButton(True, "Percentile :", True, str(90.0), False)

        self.histoButton.setMaximumWidth(250)
        self.minBin.setMaximumWidth(250)
        self.maxBin.setMaximumWidth(250)
        self.nbBin.setMaximumWidth(250)
        self.percentile.setMaximumWidth(250)

        self.histoButton.clicked.connect(self.histo_Function)

        self.buttonLayout.addWidget(self.minBin)
        self.buttonLayout.addWidget(self.maxBin)
        self.buttonLayout.addWidget(self.nbBin)
        self.buttonLayout.addWidget(self.percentile)
        self.buttonLayout.addWidget(self.histoButton)
        self.setLayout(self.mainLayout)

    def histo_Function(self):
        array = self.Data_list[self.imageSelection.currentRow()]
        minBin = float(str(self.minBin.lineEdit.text()))
        maxBin = float(str(self.maxBin.lineEdit.text()))
        nbBin = int(str(self.nbBin.lineEdit.text()))

        arrayM = np.ma.masked_outside(array, minBin + (2 * (maxBin - minBin) / (nbBin)),
                                      maxBin - (2 * (maxBin - minBin) / (nbBin)), copy=True)
        mean_value = stat.mstats.tmean(array, (
        minBin + (2 * (maxBin - minBin) / (nbBin)), maxBin - (2 * (maxBin - minBin) / (nbBin))),
                                       inclusive=(False, True))
        mode_value = stat.mstats.mode(arrayM, axis=None)
        variation_coeff = stat.mstats.variation(arrayM, axis=None)

        histo = IP.histogram(array, minBin, maxBin, nbBin)

        text = 'Mean Value : ' + str(mean_value) + '\n'
        text += 'Mode Position : ' + str(mode_value) + '\n'
        text += 'Variation Coeff : ' + str(variation_coeff) + '\n'
        text += '________Histogram________' + '\n'
        text += 'Occurance ||    Bins        ' + '\n'
        for i in range(1, len(histo[0])):
            text += str(histo[0][i]) + '\t' + str(histo[1][i]) + '\n'

        self.hide_all_button()
        self.txtDisplay = qt.QTextEdit("Results")
        self.txtDisplay.setText(text)
        self.txtDisplay.setMaximumWidth(250)
        self.buttonLayout.addWidget(self.txtDisplay)

        self.plt.show()
        self.plt.clear()
        pen_act = pg.mkPen((0, 200, 0, 255), width=1)
        inf1 = pg.InfiniteLine(movable=True, pos=str(mean_value), angle=90)
        # inf2 = pg.InfiniteLine(movable=True,pos= str(mode_value), angle=90)
        self.p1 = self.plt.addPlot(row=0, col=0, title="Histogram",
                                   labels={'left': "Voxel Number ", 'bottom': "Voxels Value"})
        self.p1.addItem(inf1, title='Mean Value')
        # self.p1.addItem(inf2,title='Mode Value')
        self.curve1 = self.p1.plot(histo[1][1:-1], histo[0][1:], pen=pen_act, stepMode=False, fillLevel=0)
        self.plt.repaint()

    def filter_anidiff_GUI(self):

        self.hide_all_button()

        self.filterButton = qt.QPushButton("Filter")
        self.time_step = LabelEditAndButton(True, "Time Step : ", True, str(0.06), False)
        self.conductance = LabelEditAndButton(True, "Conductance : ", True, str(9.0), False)
        self.iterLabel = qt.QLabel("# iterations")
        self.nbIter = SliderAndLabel(self)
        self.nbIter._setRange(1, 100)

        self.filterButton.setMaximumWidth(250)
        self.time_step.setMaximumWidth(250)
        self.conductance.setMaximumWidth(250)
        self.iterLabel.setMaximumWidth(250)
        self.nbIter.setMaximumWidth(250)

        self.filterButton.clicked.connect(self.filter_anidiff_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.time_step)
        self.buttonLayout.addWidget(self.conductance)
        self.buttonLayout.addWidget(self.iterLabel)
        self.buttonLayout.addWidget(self.nbIter)
        self.buttonLayout.addWidget(self.filterButton)
        self.setLayout(self.mainLayout)

    def filter_anidiff_Function(self, nbIter=-1):

        inputDataToFilter = self.Data_list[self.imageSelection.currentRow()]

        time_step = float(str(self.time_step.lineEdit.text()))
        conductance = float(str(self.conductance.lineEdit.text()))

        if nbIter == -1:
            nbIter = self.nbIter.value()

        ImageOut = IP.anisotropic_diffusion(np.copy(inputDataToFilter), time_step, conductance, nbIter)
        self.addImage("F_AniDif" + str(time_step) + '-' + str(conductance), ImageOut)
        self.filter_anidiff_GUI()

    def filter_recursiveGauss_GUI(self):

        self.hide_all_button()

        self.filterButton = qt.QPushButton("Filter")
        self.sigma = LabelEditAndButton(True, "Sigma : ", True, str(1.0), False)

        self.sigma.setMaximumWidth(250)

        self.filterButton.clicked.connect(self.filter_recursiveGauss_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.sigma)
        self.buttonLayout.addWidget(self.filterButton)
        self.setLayout(self.mainLayout)

    def filter_recursiveGauss_Function(self):

        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

        sigma = float(str(self.sigma.lineEdit.text()))

        ImageOut = IP.recursiveGauss(np.copy(inputDataToSeg), sigma)
        self.addImage("F_RGauss" + str(sigma), ImageOut)
        self.filter_recursiveGauss_GUI()

    def filter_MIP_GUI(self):
        self.hide_all_button()

        self.mip = qt.QComboBox()
        self.mip.addItems(['Min', 'Max'])
        self.direction = qt.QComboBox()
        self.direction.addItems(['x', 'y', 'z'])

        self.mip_thickness = LabelEditAndButton(True, "Thickness :  ", True, str(10), False)

        self.mipButton = qt.QPushButton("GO !")

        self.mipButton.clicked.connect(self.filter_mip)

        self.mip.setMaximumWidth(250)
        self.direction.setMaximumWidth(250)
        self.mip_thickness.setMaximumWidth(250)
        self.mipButton.setMaximumWidth(250)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.mip)
        self.buttonLayout.addWidget(self.direction)
        self.buttonLayout.addWidget(self.mip_thickness)
        self.buttonLayout.addWidget(self.mipButton)
        self.setLayout(self.mainLayout)

    def filter_mip(self):
        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]
        mip = self.mip.currentText()
        direction = self.direction.currentText()
        thickness = int(str(self.mip_thickness.lineEdit.text()))
        Image_out = IP.mip(inputDataToSeg, mip, direction, thickness)
        self.addImage("MIP_" + direction + '_' + mip, Image_out)

    def filter_WL_GUI(self):
        self.hide_all_button()
        self.waveletList = qt.QComboBox()
        self.waveletList.addItems(
            ['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8', 'bior3.1', 'bior3.3',
             'bior3.5', 'bior3.7', 'bior3.9', 'bior4.4', 'bior5.5', 'bior6.8', 'coif1', 'coif2', 'coif3', 'coif4',
             'coif5', 'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db11', 'db12', 'db13',
             'db14', 'db15', 'db16', 'db17', 'db18', 'db19', 'db20', 'dmey', 'haar', 'rbio1.1', 'rbio1.3', 'rbio1.5',
             'rbio2.2', 'rbio2.4', 'rbio2.6', 'rbio2.8', 'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7', 'rbio3.9',
             'rbio4.4', 'rbio5.5', 'rbio6.8', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'sym9', 'sym10',
             'sym11', 'sym12', 'sym13', 'sym14', 'sym15', 'sym16', 'sym17', 'sym18', 'sym19', 'sym20'])

        self.filterButton = qt.QPushButton("Filter")
        self.levels = LabelEditAndButton(True, "Levels : ", True, str(2), False)
        self.alpha = LabelEditAndButton(True, "Alpha : ", True, str(2), False)

        self.alpha.setMaximumWidth(250)
        self.levels.setMaximumWidth(250)

        self.filterButton.clicked.connect(self.filter_WL)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.waveletList)
        self.buttonLayout.addWidget(self.levels)
        self.buttonLayout.addWidget(self.alpha)
        self.buttonLayout.addWidget(self.filterButton)
        self.setLayout(self.mainLayout)

    def filter_WL(self):

        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

        wavelet = self.waveletList.currentText()

        levels = int(str(self.levels.lineEdit.text()))
        alpha = float(str(self.alpha.lineEdit.text()))

        Image_Out = IP.dwt_denoise(np.copy(inputDataToSeg), wavelet, levels, alpha)
        self.addImage("F_Wavelet_" + wavelet + '_' + str(alpha) + '_' + str(levels), Image_Out)
        self.filter_WL_GUI()

    def filter_median_GUI(self):

        self.hide_all_button()

        self.filterButton = qt.QPushButton("Filter")
        self.radius = LabelEditAndButton(True, "Radius : ", True, str(2.0), False)

        self.flag2d = qt.QComboBox()
        self.flag2d.addItem("2D")
        self.flag2d.addItem("3D")

        self.radius.setMaximumWidth(250)

        self.filterButton.clicked.connect(self.filter_median_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.radius)
        self.buttonLayout.addWidget(self.filterButton)
        self.buttonLayout.addWidget(self.flag2d)
        self.setLayout(self.mainLayout)

    def filter_median_Function(self):

        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]
        radius = float(str(self.radius.lineEdit.text()))

        if self.flag2d.currentIndex() == 0:
            flag2d = True
        else:
            flag2d = False

        ImageOut = IP.median(np.copy(inputDataToSeg), int(radius), flag2d)
        self.addImage("F_median" + str(radius), ImageOut)
        self.filter_median_GUI()

    def zero_Crossing_GUI(self):

        self.hide_all_button()

        self.edgeButton = qt.QPushButton("Edge")
        self.varGauss = LabelEditAndButton(True, "Variance Gaussian:", True, str(3.0), False)
        self.maxErrorGauss = LabelEditAndButton(True, "Max Error Gaussian: ", True, str(0.01), False)

        self.edgeButton.clicked.connect(self.zero_Crossing_Function)

        self.varGauss.setMaximumWidth(250)
        self.maxErrorGauss.setMaximumWidth(250)
        self.edgeButton.setMaximumWidth(250)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.varGauss)
        self.buttonLayout.addWidget(self.maxErrorGauss)
        self.buttonLayout.addWidget(self.edgeButton)

        self.setLayout(self.mainLayout)

    def zero_Crossing_Function(self):

        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

        var = float(str(self.varGauss.lineEdit.text()))
        maxErr = float(str(self.maxErrorGauss.lineEdit.text()))
        ImageOut = IP.zero_crossing(np.copy(inputDataToSeg), var, maxErr)

        self.addImage("E_Zero_" + str(var) + "_" + str(maxErr), ImageOut)

        self.zero_Crossing_GUI()

    def edge_Canny_GUI(self):

        self.hide_all_button()
        self.edgeButton = qt.QPushButton("Edge")

        self.varGauss = LabelEditAndButton(True, "Variance Gaussian:", True, str(3.0), False)
        self.maxErrorGauss = LabelEditAndButton(True, "Max Error Gaussian: ", True, str(0.01), False)

        self.varGauss.setMaximumWidth(250)
        self.maxErrorGauss.setMaximumWidth(250)
        self.edgeButton.setMaximumWidth(250)

        self.edgeButton.clicked.connect(self.edge_Canny_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.varGauss)
        self.buttonLayout.addWidget(self.maxErrorGauss)

        self.buttonLayout.addWidget(self.edgeButton)
        self.setLayout(self.mainLayout)

    def edge_Canny_Function(self):
        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]
        var = float(str(self.varGauss.lineEdit.text()))
        maxErr = float(str(self.maxErrorGauss.lineEdit.text()))
        minT = 0.0
        maxT = 0.0
        ImageOut = IP.cannyEdge(np.copy(inputDataToSeg), var, maxErr, minT, maxT)
        self.addImage("E_Canny_" + str(var) + "_" + str(maxErr), ImageOut)
        self.edge_Canny_GUI()

    def edge_GradGauss_GUI(self):

        self.hide_all_button()

        self.filterButton = qt.QPushButton("Filter")
        self.sigma = LabelEditAndButton(True, "Sigma : ", True, str(1.0), False)

        self.sigma.setMaximumWidth(250)

        self.filterButton.clicked.connect(self.filter_GradGauss_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.sigma)
        self.buttonLayout.addWidget(self.filterButton)
        self.setLayout(self.mainLayout)

    def filter_GradGauss_Function(self):

        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

        if (inputDataToSeg.shape[0] > 4) and (inputDataToSeg.shape[1] > 4) and (inputDataToSeg.shape[2] > 4):

            sigma = float(str(self.sigma.lineEdit.text()))

            ImageOut = IP.gradGauss(np.copy(inputDataToSeg), sigma)
            self.addImage("F_GradGauss" + str(sigma), ImageOut)
        else:
            msgBox = qt.QMessageBox(self)
            msgBox.setText('Image Dimension need to be larger than 4 x 4 x 4 ')
            msgBox.exec_()

        self.filter_GradGauss_GUI()

    def filter_Sigmo_GUI(self):

        self.hide_all_button()

        self.filterButton = qt.QPushButton("Filter")
        self.alpha = LabelEditAndButton(True, "Alpha : ", True, str(0.1), False)
        self.beta = LabelEditAndButton(True, "Beta: ", True, str(0.3), False)

        self.alpha.setMaximumWidth(250)
        self.beta.setMaximumWidth(250)

        self.filterButton.clicked.connect(self.filter_Sigmo_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.alpha)
        self.buttonLayout.addWidget(self.beta)
        self.buttonLayout.addWidget(self.filterButton)
        self.setLayout(self.mainLayout)

    def filter_Sigmo_Function(self):
        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

        alpha = float(str(self.alpha.lineEdit.text()))
        beta = float(str(self.beta.lineEdit.text()))

        ImageOut = IP.sigmo(np.copy(inputDataToSeg), alpha, beta)
        self.addImage("F_Sigmo" + str(alpha) + '-' + str(beta), ImageOut)

        self.filter_Sigmo_GUI()

    """ Segmentation  """

    def _rgc_segmentation_GUI(self):
        self.hide_all_button()

        self.segmentButtonrg = qt.QPushButton("Segment")
        self.rg_radius = LabelEditAndButton(True, "Radius Of Confidence: ", True, str(1), False)
        self.rg_mul = LabelEditAndButton(True, "Standart Deviation Multiplier : ", True, str(0.1), False)
        self.rg_Iter = LabelEditAndButton(True, "Number of Iteration : ", True, str(1), False)

        self.flag_indep = qt.QCheckBox("Segment Each Seed has Independent:")
        self.flag_indep.setCheckState(True)

        self.segmentButtonrg.setMaximumWidth(250)
        self.rg_radius.setMaximumWidth(250)
        self.rg_mul.setMaximumWidth(250)
        self.rg_Iter.setMaximumWidth(250)
        self.flag_indep.setMaximumWidth(250)
        self.segmentButtonrg.clicked.connect(self.segment_rgc_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.flag_indep)
        self.buttonLayout.addWidget(self.rg_radius)
        self.buttonLayout.addWidget(self.rg_mul)
        self.buttonLayout.addWidget(self.rg_Iter)
        self.buttonLayout.addWidget(self.segmentButtonrg)
        self.setLayout(self.mainLayout)

    def _rg_segmentation_GUI(self):

        self.hide_all_button()

        self.segmentButton = qt.QPushButton("Segment")
        self.rg_tol_min = LabelEditAndButton(True, "Threathold Min : ", True, str(1.0), False)
        self.rg_tol_max = LabelEditAndButton(True, "Threathold Max : ", True, str(2.0), False)

        self.segmentButton.setMaximumWidth(250)
        self.rg_tol_min.setMaximumWidth(250)
        self.rg_tol_max.setMaximumWidth(250)

        self.segmentButton.clicked.connect(self.segment_rg_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.rg_tol_min)
        self.buttonLayout.addWidget(self.rg_tol_max)
        self.buttonLayout.addWidget(self.segmentButton)
        self.setLayout(self.mainLayout)

    def segment_rgc_Function(self, seedListToSegment=-1):

        if seedListToSegment == -1:
            seedListToSegment = []

            for direction in self.ItemsLists[self.imageSelection.currentRow()]['Seeds']:
                for seed in self.ItemsLists[self.imageSelection.currentRow()]['Seeds'][direction]:
                    if seed not in seedListToSegment:
                        seedListToSegment.append(seed)

        if (len(seedListToSegment) > 0):
            inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

            radius = int(str(self.rg_radius.lineEdit.text()))
            multi = float(str(self.rg_mul.lineEdit.text()))
            iterN = int(str(self.rg_Iter.lineEdit.text()))
            if not bool(self.flag_indep.checkState()):
                ImageOut = IP.SegConnectedThresholdC(np.copy(inputDataToSeg), radius, multi, iterN, seedListToSegment)
            else:
                ImageOut = IP.SegConnectedThresholdC(np.copy(inputDataToSeg), radius, multi, iterN,
                                                     [seedListToSegment[0]])
                for seed in seedListToSegment[1:]:
                    ImageOut += IP.SegConnectedThresholdC(np.copy(inputDataToSeg), radius, multi, iterN, [seed])

            self.addImage("Seg " + str(radius) + '-' + str(multi), ImageOut)

        else:
            msgBox = qt.QMessageBox(self)
            msgBox.setText('Please Select a starting Point ')
            msgBox.exec_()

        self._rgc_segmentation_GUI()

    def _wh_segmentation_GUI(self):

        self.hide_all_button()
        self.segmentButton = qt.QPushButton("Segment")
        self.wh_level = LabelEditAndButton(True, "Level : ", True, str(1), False)
        self.labelLine = qt.QLabel("Watershed Line")
        self.markLine = qt.QCheckBox()
        self.connected = qt.QLabel("Fully Connected")
        self.connectedCb = qt.QCheckBox()

        self.segmentButton.setMaximumWidth(250)
        self.wh_level.setMaximumWidth(250)
        self.labelLine.setMaximumWidth(250)
        self.markLine.setMaximumWidth(250)
        self.connected.setMaximumWidth(250)
        self.connectedCb.setMaximumWidth(250)

        self.segmentButton.clicked.connect(self.segment_wh_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.wh_level)
        self.buttonLayout.addWidget(self.labelLine)
        self.buttonLayout.addWidget(self.markLine)
        self.buttonLayout.addWidget(self.connected)
        self.buttonLayout.addWidget(self.connectedCb)
        self.buttonLayout.addWidget(self.segmentButton)
        self.setLayout(self.mainLayout)

    def segment_wh_Function(self, level=-1):
        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]
        if level == -1:
            level = float(self.wh_level.lineEdit.text())

        markLine = bool(self.connectedCb.checkState())
        if markLine == 2:
            markLine = True
        else:
            markLine = False

        connected = bool(self.connectedCb.checkState())

        if connected == 2:
            connected = True
        else:
            connected = False

        OutputData = IP.SegWatershed(inputDataToSeg, level, markLine, connected)
        self.addImage("Water_" + str(level) + '_' + self.Name_list[self.imageSelection.currentRow()], OutputData)
        self._wh_segmentation_GUI()

    def segment_rg_Function(self, minTh=-1, maxTh=-1, seedListToSegment=-1):

        if seedListToSegment == - 1:
            seedListToSegment = []

            for direction in self.ItemsLists[self.imageSelection.currentRow()]['Seeds']:
                for seed in self.ItemsLists[self.imageSelection.currentRow()]['Seeds'][direction]:
                    seedListToSegment.append(seed)

        if (len(seedListToSegment) > 0):
            inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]
            if (minTh == False) and (maxTh == -1):
                minTh = float(str(self.rg_tol_min.lineEdit.text()))
                maxTh = float(str(self.rg_tol_max.lineEdit.text()))

            ImageOut = IP.SegConnectedThreshold(np.copy(inputDataToSeg), minTh, maxTh, seedListToSegment)
            self.addImage("Seg " + str(minTh) + '-' + str(maxTh), ImageOut)
            return ImageOut
        else:
            msgBox = qt.QMessageBox(self)
            msgBox.setText('Please Select a starting Point ')
            msgBox.exec_()

        self._rg_segmentation_GUI()

    def _fm_segmentation_GUI(self):

        self.hide_all_button()

        self.segmentButton = qt.QPushButton("Segment")

        self.stopValue = LabelEditAndButton(True, "Stop Value: ", True, str(1.0), False)

        self.segmentButton.setMaximumWidth(250)
        self.stopValue.setMaximumWidth(250)

        self.segmentButton.clicked.connect(self.segment_fm_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.stopValue)
        self.buttonLayout.addWidget(self.segmentButton)
        self.setLayout(self.mainLayout)

    def segment_fm_Function(self):

        seedListToSegment = []

        directions = self.ItemsList[self.imageSelection.currentRow()]['Seeds']

        for direction in directions:

            for seed in direction:
                if len(seedListToSegment) == 0:
                    seedListToSegment = seed
                else:
                    seedListToSegment.append(seed)

        if (len(seedListToSegment) > 0):

            inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

            stopValue = float(str(self.stopValue.lineEdit.text()))
            ImageOut = IP.FastMarching(np.copy(inputDataToSeg), stopValue, seedListToSegment)
            self.addImage("FM_ " + str(stopValue), ImageOut)

        else:
            msgBox = qt.QMessageBox(self)
            msgBox.setText('Please Select a starting Point ')
            msgBox.exec_()

        self._fm_segmentation_GUI()

    def _th_segmentation_GUI(self):
        self.hide_all_button()

        self.segmentButton = qt.QPushButton("Segment")
        self.min_th = LabelEditAndButton(True, "Minimum Threshold: ", True, str(0.0), False)
        self.max_th = LabelEditAndButton(True, "Maximum Threshold: ", True, str(1.0), False)

        self.segmentButton.setMaximumWidth(250)
        self.min_th.setMaximumWidth(250)
        self.max_th.setMaximumWidth(250)

        self.segmentButton.clicked.connect(self.segment_th_Function)

        self.buttonLayout.setAlignment(qt.Qt.AlignTop)
        self.buttonLayout.addWidget(self.min_th)
        self.buttonLayout.addWidget(self.max_th)
        self.buttonLayout.addWidget(self.segmentButton)
        self.setLayout(self.mainLayout)

    def segment_th_Function(self):

        inputDataToSeg = self.Data_list[self.imageSelection.currentRow()]

        min_th = float(str(self.min_th.lineEdit.text()))
        max_th = float(str(self.max_th.lineEdit.text()))
        ImageOut = IP.Threshold(np.copy(inputDataToSeg), min_th, max_th)
        self.addImage("Th_ " + str(min_th) + str(max_th), ImageOut)

        self._th_segmentation_GUI()

    def _sd_segmentation_GUI(self):

        self.hide_all_button()

        self.segmentButton = qt.QPushButton("Segment")
        self.label1 = qt.QLabel("Initial level set Image")
        self.InitialLS = qt.QComboBox()

        self.radius = LabelEditAndButton(False, "Radius", True, "3.0", False)
        self.label2 = qt.QLabel("Edge Potential Map")
        self.EdgePotMap = qt.QComboBox()
        self.maxRMSError = LabelEditAndButton(True, "Max RMS Error: ", True, str(1.0), False)
        self.propaScaling = LabelEditAndButton(True, "Propagation Scaling: ", True, str(1.0), False)
        self.curvScaling = LabelEditAndButton(True, "Curvature Scaling: ", True, str(1.0), False)

        self.iterLabel = qt.QLabel("# iterations")
        self.nbIter = SliderAndLabel(self)
        self.nbIter._setRange(1, 100000)

        self.segmentButton.setMaximumWidth(250)
        self.label1.setMaximumWidth(250)
        self.InitialLS.setMaximumWidth(250)
        self.label2.setMaximumWidth(250)
        self.EdgePotMap.setMaximumWidth(250)
        self.maxRMSError.setMaximumWidth(250)
        self.propaScaling.setMaximumWidth(250)
        self.curvScaling.setMaximumWidth(250)
        self.nbIter.setMaximumWidth(250)
        self.radius.setMaximumWidth(250)

        self.InitialLS.addItem("Starting Seed Points")
        self.InitialLS.addItems(self.Name_list)
        self.EdgePotMap.addItems(self.Name_list)

        self.segmentButton.clicked.connect(self.segment_sd_Function)
        self.InitialLS.currentIndexChanged.connect(self._radius)

        self.buttonLayout.addWidget(self.label1)
        self.buttonLayout.addWidget(self.InitialLS)
        self.buttonLayout.addWidget(self.radius)
        self.buttonLayout.addWidget(self.label2)
        self.buttonLayout.addWidget(self.EdgePotMap)
        self.buttonLayout.addWidget(self.maxRMSError)
        self.buttonLayout.addWidget(self.propaScaling)
        self.buttonLayout.addWidget(self.curvScaling)
        self.buttonLayout.addWidget(self.iterLabel)
        self.buttonLayout.addWidget(self.nbIter)
        self.buttonLayout.addWidget(self.segmentButton)


    def _snake_chan_vese_GUI(self):

        self.hide_all_button()

        self.segmentButton = qt.QPushButton("Segment")
        self.iter = LabelEditAndButton(False, "Iteration", True, "20", False)



        self.label1 = qt.QLabel("Initial level set Image")
        self.InitialLS = qt.QComboBox()

        self.radius = LabelEditAndButton(False, "Radius", True, "3.0", False)
        self.label2 = qt.QLabel("Edge Potential Map")
        self.EdgePotMap = qt.QComboBox()
        self.maxRMSError = LabelEditAndButton(True, "Max RMS Error: ", True, str(1.0), False)
        self.propaScaling = LabelEditAndButton(True, "Propagation Scaling: ", True, str(1.0), False)
        self.curvScaling = LabelEditAndButton(True, "Curvature Scaling: ", True, str(1.0), False)



    def segment_sd_Function(self):

        maxRMSError = float(str(self.maxRMSError.lineEdit.text()))
        propaScaling = float(str(self.propaScaling.lineEdit.text()))
        curvScaling = float(str(self.curvScaling.lineEdit.text()))
        nbIter = self.nbIter.value()

        if self.InitialLS.currentIndex() == 0:

            seedListToSegment = []

            directions = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']

            for direction in directions:
                for seed in self.ItemsLists[self.imageSelection.currentRow()]['Seeds'][direction]:
                    seedListToSegment.append(seed)

            if (len(seedListToSegment) > 0):

                radius = float(str(self.radius.lineEdit.text()))
                if radius <= 0.0:
                    radius = 1

                EdgePotMap = self.Data_list[self.EdgePotMap.currentIndex()]
                ImageOut = IP.ShapeDetectionLS([], EdgePotMap, maxRMSError, propaScaling, curvScaling, nbIter, radius,
                                               seedListToSegment)
            else:

                msgBox = qt.QMessageBox(self)
                msgBox.setText('Please Select a starting Point ')
                msgBox.exec_()
        else:

            InitLS = self.Data_list[self.InitialLS.currentIndex() - 1]
            EdgePotMap = self.Data_list[self.EdgePotMap.currentIndex()]
            ImageOut = IP.ShapeDetectionLS(InitLS, EdgePotMap, maxRMSError, propaScaling, curvScaling, nbIter, 0, [])

        self.addImage("SD_ " + str(propaScaling) + "-" + str(curvScaling) + "-" + str(nbIter), ImageOut)
        self._sd_segmentation_GUI()

    def _geo_segmentation_GUI(self):

        self.hide_all_button()

        self.segmentButton = qt.QPushButton("Segment")
        self.label1 = qt.QLabel("Initial level set Image")
        self.InitialLS = qt.QComboBox()

        self.radius = LabelEditAndButton(False, "Radius", True, "3.0", False)
        self.label2 = qt.QLabel("Edge Potential Map")
        self.EdgePotMap = qt.QComboBox()
        self.maxRMSError = LabelEditAndButton(True, "Max RMS Error: ", True, str(1.0), False)
        self.AdvectionScaling = LabelEditAndButton(True, "Advection Scaling: ", True, str(1.0), False)
        self.propaScaling = LabelEditAndButton(True, "Propagation Scaling: ", True, str(1.0), False)
        self.curvScaling = LabelEditAndButton(True, "Curvature Scaling: ", True, str(1.0), False)

        self.iterLabel = qt.QLabel("# iterations")
        self.nbIter = SliderAndLabel(self)
        self.nbIter._setRange(1, 100000)

        self.AdvectionScaling.setMaximumWidth(250)
        self.segmentButton.setMaximumWidth(250)
        self.label1.setMaximumWidth(250)
        self.InitialLS.setMaximumWidth(250)
        self.label2.setMaximumWidth(250)
        self.EdgePotMap.setMaximumWidth(250)
        self.maxRMSError.setMaximumWidth(250)
        self.propaScaling.setMaximumWidth(250)
        self.curvScaling.setMaximumWidth(250)
        self.nbIter.setMaximumWidth(250)
        self.radius.setMaximumWidth(250)

        self.InitialLS.addItem("Starting Seed Points")
        self.InitialLS.addItems(self.Name_list)
        self.EdgePotMap.addItems(self.Name_list)

        self.segmentButton.clicked.connect(self.segment_geo_Function)
        self.InitialLS.currentIndexChanged.connect(self._radius)

        self.buttonLayout.addWidget(self.label1)
        self.buttonLayout.addWidget(self.InitialLS)
        self.buttonLayout.addWidget(self.radius)
        self.buttonLayout.addWidget(self.label2)
        self.buttonLayout.addWidget(self.EdgePotMap)
        self.buttonLayout.addWidget(self.maxRMSError)
        self.buttonLayout.addWidget(self.AdvectionScaling)
        self.buttonLayout.addWidget(self.propaScaling)
        self.buttonLayout.addWidget(self.curvScaling)
        self.buttonLayout.addWidget(self.iterLabel)
        self.buttonLayout.addWidget(self.nbIter)
        self.buttonLayout.addWidget(self.segmentButton)

    def segment_geo_Function(self):

        maxRMSError = float(str(self.maxRMSError.lineEdit.text()))
        propaScaling = float(str(self.propaScaling.lineEdit.text()))
        curvScaling = float(str(self.curvScaling.lineEdit.text()))
        advScaling = float(str(self.AdvectionScaling.lineEdit.text()))

        nbIter = self.nbIter.value()

        if self.InitialLS.currentIndex() == 0:

            seedListToSegment = []

            directions = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']

            for direction in directions:

                for seed in self.ItemsLists[self.imageSelection.currentRow()]['Seeds'][direction]:
                    seedListToSegment.append(seed)

            if (len(seedListToSegment) > 0):

                radius = float(str(self.radius.lineEdit.text()))
                if radius <= 0.0:
                    radius = 1

                EdgePotMap = self.Data_list[self.EdgePotMap.currentIndex()]
                ImageOut = IP.GeoDetectionLS([], EdgePotMap, maxRMSError, propaScaling, curvScaling, advScaling, nbIter,
                                             radius, seedListToSegment)
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Please Select a starting Point ')
                msgBox.exec_()
        else:

            InitLS = self.Data_list[self.InitialLS.currentIndex() - 1]
            EdgePotMap = self.Data_list[self.EdgePotMap.currentIndex()]
            ImageOut = IP.GeoDetectionLS(InitLS, EdgePotMap, maxRMSError, propaScaling, curvScaling, advScaling, nbIter,
                                         0, [])

        self.addImage("GEO_ " + str(propaScaling), ImageOut)
        self._geo_segmentation_GUI()

    """ Other  """

    def _math_GUI(self):

        self.hide_all_button()

        self.Image1 = qt.QComboBox()
        self.Image2 = qt.QComboBox()
        self.Operator = qt.QComboBox()
        self.constante1 = LabelEditAndButton(False, "", True, "1", False)
        self.constante2 = LabelEditAndButton(False, "", True, "1", False)
        self.mathButton = qt.QPushButton("Compute")

        self.Image1.addItem(".")
        self.Image1.addItem("Constant")
        self.Image1.addItems(self.Name_list)

        self.Image2.addItem(".")
        self.Image2.addItem("Constant")
        self.Image2.addItems(self.Name_list)

        self.Operator.addItem("+")
        self.Operator.addItem("-")
        self.Operator.addItem("x")
        self.Operator.addItem("/")
        self.Operator.addItem("^")
        self.Operator.addItem("e")
        self.Operator.addItem("log")
        self.Operator.addItem("log10")
        self.Operator.addItem("&")

        self.constante1.setFixedSize(100, 40)

        self.constante2.setFixedSize(100, 40)

        self.Image1.currentIndexChanged.connect(self._constante1)
        self.Image2.currentIndexChanged.connect(self._constante2)
        self.mathButton.clicked.connect(self.mathFunction)

        self.buttonLayout.addWidget(self.Image1)
        self.buttonLayout.addWidget(self.constante1)
        self.buttonLayout.addWidget(self.Operator)
        self.buttonLayout.addWidget(self.Image2)
        self.buttonLayout.addWidget(self.constante2)
        self.buttonLayout.addWidget(self.mathButton)

        self.constante1.hide()
        self.constante2.hide()

    def _segFromContour(self):

        if len(self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction0"]) != 0:
            self._smoothContours()
            self._smoothContours()

            x_size = self.Data_list[self.imageSelection.currentRow()].shape[0]
            y_size = self.Data_list[self.imageSelection.currentRow()].shape[1]
            z_size = self.Data_list[self.imageSelection.currentRow()].shape[2]

            mask_Out = np.zeros((x_size, y_size, z_size))

            new_arr = []
            for i, curve in enumerate(self.ItemsLists[self.imageSelection.currentRow()]["Poly"]["Direction0"]):
                if len(curve) != 0:
                    for coord in curve:
                        if len(coord) != 0:
                            mask_Out[:, coord[1], coord[0]] = 1
                            new_arr.append(coord)

            mask_Out = IP.morpho('Fill', mask_Out, 0)

            self.addImage('SegContour', mask_Out * self.Data_list[self.imageSelection.currentRow()])

    def mathFunction(self):

        FlagOut = 0
        NameOperator = self.Operator.currentText()
        if self.Image1.currentText() == 'Constant':
            Operand1 = float(self.constante1.lineEdit.text())
            Name_Oper1 = self.constante1.lineEdit.text() + ' '
        elif self.Image1.currentText() == '.':
            Operand1 = None
            Name_Oper1 = '.'
        else:
            Operand1 = self.Data_list[self.Image1.currentIndex() - 2]
            Name_Oper1 = 'Im '

        if self.Image2.currentText() == 'Constant':
            Operand2 = float(self.constante2.lineEdit.text())
            Name_Oper2 = ' ' + self.constante2.lineEdit.text()

        elif self.Image2.currentText() == '.':
            Operand2 = None
            Name_Oper2 = '.'
        else:
            Operand2 = self.Data_list[self.Image2.currentIndex() - 2]
            Name_Oper2 = ' Im'

        if self.Operator.currentText() == '+':
            if (Operand1 != None and Operand2 != None) and (
                    self.Image1.currentText() != 'Constant' or self.Image2.currentText() != 'Constant'):
                ImageOut = Operand1 + Operand2
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for addition ')
                msgBox.exec_()

        if self.Operator.currentText() == '-':
            ImageOut = Operand1 - Operand2
            FlagOut = 1

        if self.Operator.currentText() == '/':
            if (Operand1 != None and Operand2 != None) and (
                    self.Image1.currentText() != 'Constant' or self.Image2.currentText() != 'Constant'):
                ImageOut = np.divide(Operand1, Operand2)
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for division ')
                msgBox.exec_()

        if self.Operator.currentText() == 'x':
            if (Operand1 != None and Operand2 != None) and (
                    self.Image1.currentText() != 'Constant' or self.Image2.currentText() != 'Constant'):
                ImageOut = np.multiply(Operand1, Operand2)
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for multiplication ')
                msgBox.exec_()

        if self.Operator.currentText() == '^':
            if (Operand1 != None and Operand2 != None) and (
                    self.Image1.currentText() != 'Constant' or self.Image2.currentText() != 'Constant'):
                ImageOut = np.power(Operand1, Operand2)
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for power ')
                msgBox.exec_()

        if self.Operator.currentText() == 'e':
            if (Operand2 != None) and (self.Image2.currentText() != 'Constant'):
                ImageOut = np.exp(Operand2)
                FlagOut = 1
            else:
                msgBox.setText('Missing operand for expo ')
                msgBox.exec_()

        if self.Operator.currentText() == 'log':
            if (Operand2 != None) and (self.Image2.currentText() != 'Constant'):
                ImageOut = np.log(Operand2)
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for log ')
                msgBox.exec_()

        if self.Operator.currentText() == 'log10':
            if (Operand2 != None) and (self.Image2.currentText() != 'Constant'):
                ImageOut = np.log10(Operand2)
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for log10 ')
                msgBox.exec_()

        if self.Operator.currentText() == '&':
            if (Operand1 != None and Operand2 != None) and (
                    self.Image1.currentText() != 'Constant' or self.Image2.currentText() != 'Constant'):
                ImageOut = Operand1 + Operand2
                ImageOut[ImageOut == 2] = 1
                FlagOut = 1
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Missing operand for +M')
                msgBox.exec_()

        if FlagOut == 1:
            self.addImage(Name_Oper1 + NameOperator + Name_Oper2, ImageOut)

        self._math_GUI()

    def _equalize(self):

        zones = []
        directions = self.ItemsLists[self.imageSelection.currentRow()]['Zones']
        for direction in directions:

            for zone in self.ItemsLists[self.imageSelection.currentRow()]['Zones'][direction]:
                zones.append(zone)

        if (len(zones) == 0):
            msgBox = qt.QMessageBox(self)
            msgBox.setText('Please Select One or two uniform zone ')
            msgBox.exec_()
        else:

            zones_in = IP.returnZonesForEqualize(zones, self.Data_list[self.imageSelection.currentRow()])
            Image = IP.equalize(np.copy(self.Data_list[self.imageSelection.currentRow()]), len(zones_in), zones_in)
            'Eq_' + str(self.imageSelection.count())
            self.addImage('Eq_' + str(self.imageSelection.count()), Image)

    def _equalizePath(self):

        pts = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction0']
        Image = self.Data_list[self.imageSelection.currentRow()]

        x0 = pts[0][0]
        x1 = pts[1][0]

        y0 = pts[0][1]
        y1 = pts[1][1]

        z0 = pts[0][2]
        z1 = pts[1][2]

        if z0 > 0:
            x0 = pts[0][0]
            x1 = pts[0][0]

            y0 = pts[0][1]
            y1 = pts[0][1]

            z0 = pts[0][2]
            z1 = pts[1][2]

        indexPts = 1
        posXY = []
        for z in range(0, Image.shape[0]):

            x1 = 0
            y1 = 0
            index = 0

            for pt in pts:
                if pt[2] < z:
                    index += 1
                else:
                    break

            if (index < len(pts)):
                x0 = float(pts[index - 1][0])
                x1 = float(pts[index][0])
                y0 = float(pts[index - 1][1])
                y1 = float(pts[index][1])
                z0 = float(pts[index - 1][2])
                z1 = float(pts[index][2])

            else:

                x0 = float(pts[index - 1][0])
                x1 = float(pts[index - 1][0])
                y0 = float(pts[index - 1][1])
                y1 = float(pts[index - 1][1])
                z0 = float(pts[index - 1][2])
                z1 = float(pts[index - 1][2])

            if z0 == z1:
                xC = x0
                yC = y0
            else:
                bX = (x1 * z0 - x0 * z1) / (z0 - z1)
                aX = (x0 / z0) - (bX / z0)
                bY = (y1 * z0 - y0 * z1) / (z0 - z1)
                aY = (y0 / z0) - (bY / z0)

                xC = aX * z + bX
                yC = aY * z + bY

            posXY.append([yC, xC, z])

        Image = IP.equalizeAlongLine(Image, posXY)
        self.addImage('Equalize_', Image)

    def _equalizeHisto(self):

        NewVol = IP.equalizeHisto(self.Data_list[self.imageSelection.currentRow()])
        # self.addImage('EqHisto_' + str(self.imageSelection.count()),NewVol)

    def _crop(self):

        zones = []
        directions = self.ItemsLists[self.imageSelection.currentRow()]['Zones']

        for direction in directions:
            for zone in self.ItemsLists[self.imageSelection.currentRow()]['Zones'][direction]:
                zones.append(zone)

        if zones != []:
            zones_in = IP.returnZonesForEqualize(zones, self.Data_list[self.imageSelection.currentRow()])
            x1 = zones_in[-1][2]
            x2 = zones_in[-1][5]
            y1 = zones_in[-1][1]
            y2 = zones_in[-1][4]
            z1 = zones_in[-1][0]
            z2 = zones_in[-1][3]

            if x1 == x2:
                Image = self.Data_list[self.imageSelection.currentRow()][z1:z2, y1:y2, :]
                self.addImage('Cr_' + str(self.imageSelection.count()), Image)

            if y1 == y2:
                Image = self.Data_list[self.imageSelection.currentRow()][z1:z2, :, x1:x2]
                self.addImage('Cr_' + str(self.imageSelection.count()), Image)
            if z1 == z2:
                Image = self.Data_list[self.imageSelection.currentRow()][:, y1:y2, x1:x2]
                self.addImage('Cr_' + str(self.imageSelection.count()), Image)



        else:

            msgBox = qt.QMessageBox(self)
            msgBox.setText('Please Select at Least a zone or a contour')
            msgBox.exec_()

    def _fill(self):

        zones = []
        directions = self.ItemsLists[self.imageSelection.currentRow()]['Zones']

        for direction in directions:
            for zone in self.ItemsLists[self.imageSelection.currentRow()]['Zones'][direction]:
                zones.append(zone)

        if zones != []:
            zones_in = IP.returnZonesForEqualize(zones, self.Data_list[self.imageSelection.currentRow()])
            x1 = zones_in[-1][2]
            x2 = zones_in[-1][5]
            y1 = zones_in[-1][1]
            y2 = zones_in[-1][4]
            z1 = zones_in[-1][0]
            z2 = zones_in[-1][3]

            if x1 == x2:
                self.Data_list[self.imageSelection.currentRow()][z1:z2, y1:y2, :] = 1

            if y1 == y2:
                self.Data_list[self.imageSelection.currentRow()][z1:z2, y1:y2, :] = 1
            if z1 == z2:
                self.Data_list[self.imageSelection.currentRow()][z1:z2, y1:y2, :] = 1

        else:

            msgBox = qt.QMessageBox(self)
            msgBox.setText('Please Select at Least a zone or a contour')
            msgBox.exec_()

    def _contrast_GUI(self):

        self.hide_all_button()

        self.adjContrastButton = qt.QPushButton("Adjust Contrast")
        self.minIn = LabelEditAndButton(True, "Min Input", True, str(1), False)
        self.maxIn = LabelEditAndButton(True, "Max Input", True, str(1), False)
        self.minOut = LabelEditAndButton(True, "Min Output", True, str(1), False)
        self.maxOut = LabelEditAndButton(True, "Max Output", True, str(1), False)

        self.adjContrastButton.clicked.connect(self.contrast_function)
        self.minIn.setMaximumWidth(250)
        self.maxIn.setMaximumWidth(250)
        self.minOut.setMaximumWidth(250)
        self.maxOut.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.minIn)
        self.buttonLayout.addWidget(self.maxIn)
        self.buttonLayout.addWidget(self.minOut)
        self.buttonLayout.addWidget(self.maxOut)
        self.buttonLayout.addWidget(self.adjContrastButton)

    def contrast_function(self):

        image2adjust = self.Data_list[self.imageSelection.currentRow()]
        minIn = float(self.minIn.lineEdit.text())
        maxIn = float(self.maxIn.lineEdit.text())
        minOut = float(self.minOut.lineEdit.text())
        maxOut = float(self.maxOut.lineEdit.text())

        a = (maxOut - minOut) / (maxIn - minIn)
        b = maxOut - (a * maxIn)
        imageOut = a * image2adjust + b

        self.addImage('ReAdjust', imageOut)

    def _resample_GUI(self):

        self.hide_all_button()

        self.resampleButton = qt.QPushButton("Resample")
        self.InterpolatorLabel = qt.QLabel("Interpolator")
        self.Interpolator = qt.QComboBox()
        self.Interpolator.addItem("Nearest Neighbor")
        self.Interpolator.addItem("Linear")
        self.Interpolator.addItem("BSpline")
        self.Interpolator.addItem("Gaussian")
        self.Interpolator.addItem("HammingWindowedSinc")
        self.Interpolator.addItem("CosineWindowedSinc")
        self.Interpolator.addItem("WelchWindowedSinc")
        self.Interpolator.addItem("LanczosWindowedSinc")
        self.Interpolator.addItem("BlackmanWindowedSinc")

        self.Factorlabel = qt.QLabel("Resize Factor ")
        self.layoutFactor = qt.QHBoxLayout()
        self.xfactor = LabelEditAndButton(True, "Z", True, str(1), False)
        self.yfactor = LabelEditAndButton(True, "Y", True, str(1), False)
        self.zfactor = LabelEditAndButton(True, "X", True, str(1), False)

        self.resampleButton.clicked.connect(self.resample_function)

        self.resampleButton.setMaximumWidth(250)
        self.InterpolatorLabel.setMaximumWidth(250)
        self.Interpolator.setMaximumWidth(250)
        self.Factorlabel.setMaximumWidth(250)
        self.xfactor.setMaximumWidth(60)
        self.yfactor.setMaximumWidth(60)
        self.zfactor.setMaximumWidth(60)

        self.buttonLayout.addWidget(self.Interpolator)
        self.buttonLayout.addWidget(self.Factorlabel)
        self.buttonLayout.addWidget(self.zfactor)
        self.buttonLayout.addWidget(self.yfactor)
        self.buttonLayout.addWidget(self.xfactor)

        self.buttonLayout.addWidget(self.resampleButton)

    def _rotation_GUI(self):

        self.hide_all_button()

        self.resampleButton = qt.QPushButton("Resample")

        self.kRotation = LabelEditAndButton(True, "Rotation", True, str(90), False)

        self.resampleButton.clicked.connect(self.rotation_function)

        self.kRotation.setMaximumWidth(250)
        self.resampleButton.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.kRotation)
        self.buttonLayout.addWidget(self.resampleButton)

    def rotation_function(self):

        image = self.Data_list[self.imageSelection.currentRow()]

        rotation = float(self.kRotation.lineEdit.text())

        image = IP.rotation(image, rotation)
        self.addImage('Rotate', image)

    def _mozaicGUI(self):

        self.hide_all_button()

        self.mozaicButton = qt.QPushButton("Mozaic")

        self.label1 = qt.QLabel("Image 1")
        self.Im1 = qt.QComboBox()
        self.label2 = qt.QLabel("Image 2")
        self.Im2 = qt.QComboBox()

        self.Im1.addItems(self.Name_list)
        self.Im2.addItems(self.Name_list)

        self.labelDirection = qt.QLabel("Direction")
        self.direction = qt.QComboBox()

        self.direction.addItem("Vertical")
        self.direction.addItem("Horizontal")

        self.mozaicButton.clicked.connect(self.mozaic_function)

        self.mozaicButton.setMaximumWidth(250)
        self.label1.setMaximumWidth(250)
        self.Im1.setMaximumWidth(250)
        self.label2.setMaximumWidth(250)
        self.Im2.setMaximumWidth(250)
        self.labelDirection.setMaximumWidth(250)
        self.direction.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.label1)
        self.buttonLayout.addWidget(self.Im1)
        self.buttonLayout.addWidget(self.label2)
        self.buttonLayout.addWidget(self.Im2)
        self.buttonLayout.addWidget(self.labelDirection)
        self.buttonLayout.addWidget(self.direction)
        self.buttonLayout.addWidget(self.mozaicButton)

    def mozaic_function(self):

        Im1 = self.Data_list[self.Im1.currentIndex()]
        Im2 = self.Data_list[self.Im2.currentIndex()]

        if self.direction.currentIndex() == 0:
            size_z = np.max([Im1.shape[0], Im2.shape[0]])
            size_x = np.max([Im1.shape[1], Im2.shape[1]])
            size_y = Im1.shape[2] + Im2.shape[2]
            new_image = np.zeros((size_z, size_x, size_y))
            new_image[0:Im1.shape[0], 0:Im1.shape[1], 0:Im1.shape[2]] = Im1
            new_image[0:Im2.shape[0], 0:Im2.shape[1], Im1.shape[2]:Im1.shape[2] + Im2.shape[2]] = Im2

        else:
            size_z = np.max([Im1.shape[0], Im2.shape[0]])
            size_x = Im1.shape[1] + Im2.shape[1]
            size_y = np.max([Im1.shape[2], Im2.shape[2]])
            new_image = np.zeros((size_z, size_x, size_y))
            new_image[0:Im1.shape[0], 0:Im1.shape[1], 0:Im1.shape[2]] = Im1
            new_image[0:Im2.shape[0], Im1.shape[1]:Im1.shape[1] + Im2.shape[1], 0:Im2.shape[2]] = Im2

        self.addImage('Mozaic', new_image)

    def resample_function(self):

        image2resample = self.Data_list[self.imageSelection.currentRow()]
        interpolator = self.Interpolator.currentText()
        zf = float(self.zfactor.lineEdit.text())
        xf = float(self.xfactor.lineEdit.text())
        yf = float(self.yfactor.lineEdit.text())
        factor = [1.0 / zf, 1.0 / yf, 1.0 / xf]

        imageOut = IP.resizeImage(np.copy(image2resample), interpolator, factor)

        self.addImage('Resized_' + str(xf) + '_' + str(yf) + '_' + str(zf), imageOut)

    def _stack_GUI(self):

        self.hide_all_button()

        self.stackButton = qt.QPushButton("Stack")
        self.label1 = qt.QLabel("TopImage")
        self.Im1 = qt.QComboBox()
        self.index_z1_i1 = LabelEditAndButton(True, "Z1: ", True, str(0), False)
        self.index_z2_i1 = LabelEditAndButton(True, "Z2: ", True, str(10), False)
        self.label2 = qt.QLabel("Botom Image")
        self.Im2 = qt.QComboBox()
        self.index_z1_i2 = LabelEditAndButton(True, "Z1: ", True, str(0), False)
        self.index_z2_i2 = LabelEditAndButton(True, "Z2: ", True, str(10), False)

        self.stackButton.setMaximumWidth(250)
        self.label1.setMaximumWidth(250)
        self.Im1.setMaximumWidth(250)
        self.index_z1_i1.setMaximumWidth(250)
        self.index_z2_i1.setMaximumWidth(250)
        self.label2.setMaximumWidth(250)
        self.Im2.setMaximumWidth(250)
        self.index_z1_i2.setMaximumWidth(250)
        self.index_z2_i2.setMaximumWidth(250)

        self.Im1.addItems(self.Name_list)
        self.Im2.addItems(self.Name_list)

        self.stackButton.clicked.connect(self.stack_Function)

        self.buttonLayout.addWidget(self.label1)
        self.buttonLayout.addWidget(self.Im1)
        self.buttonLayout.addWidget(self.index_z1_i1)
        self.buttonLayout.addWidget(self.index_z2_i1)
        self.buttonLayout.addWidget(self.label2)
        self.buttonLayout.addWidget(self.Im2)
        self.buttonLayout.addWidget(self.index_z1_i2)
        self.buttonLayout.addWidget(self.index_z2_i2)

        self.buttonLayout.addWidget(self.stackButton)

    def _destack_GUI(self):
        self.hide_all_button()

        self.destackButton = qt.QPushButton("Destack")
        self.index_h = LabelEditAndButton(True, "Height Stack: ", True, str(60), False)
        self.index_zmin = LabelEditAndButton(True, "Z min coordinate : ", True, str(4), False)
        self.index_zmax = LabelEditAndButton(True, "Z max  coordinate: ", True, str(55), False)

        self.destackButton.setMaximumWidth(250)
        self.index_zmin.setMaximumWidth(250)
        self.index_zmax.setMaximumWidth(250)
        self.index_h.setMaximumWidth(250)

        self.destackButton.clicked.connect(self.destack_Function)

        self.buttonLayout.addWidget(self.index_h)
        self.buttonLayout.addWidget(self.index_zmin)
        self.buttonLayout.addWidget(self.index_zmax)
        self.buttonLayout.addWidget(self.destackButton)

    def _followLineGUI(self):
        self.hide_all_button()
        self.MIButton = qt.QPushButton("Plot Line Intensity")
        self.MIButton.setMaximumWidth(250)
        self.MIButton.clicked.connect(self.followLine)
        self.buttonLayout.addWidget(self.MIButton)

    def _followMaxLineGUI(self):
        self.hide_all_button()
        self.MIButton = qt.QPushButton("Plot Line max Intensity")
        self.maxDistance = LabelEditAndButton(True, "Pixel To Travel: ", True, str(1000), False)
        self.MIButton.setMaximumWidth(250)
        self.maxDistance.setMaximumWidth(250)

        self.MIButton.clicked.connect(self.followMaxLine)
        self.buttonLayout.addWidget(self.maxDistance)
        self.buttonLayout.addWidget(self.MIButton)

    def followLine(self):

        image2Folow = self.Data_list[self.imageSelection.currentRow()]
        seed1 = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction0'][0]
        seed2 = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction0'][1]
        dic = IP.folowLine(image2Folow, seed1, seed2, [0, 0, 0])

        diameter = dic['Diameter']
        distance = dic['Distance']
        Image = dic['Image']

        self.plt.show()
        self.plt.clear()
        self.plt.showGrid(x=True, y=True)
        self.plt.setLabel('left', 'Diameter')
        self.plt.setLabel('bottom', 'Distance')
        self.curve = self.plt.plot(diameter, distance)
        self.plt.repaint()
        self.addImage('FollowedLine', Image)

    def followMaxLine(self, maxDistance=-1):

        image2Folow = self.Data_list[self.imageSelection.currentRow()]
        seed1 = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction0'][0]
        seed2 = self.ItemsLists[self.imageSelection.currentRow()]['Seeds']['Direction0'][1]

        if maxDistance == -1:
            maxDistance = int(str(self.maxDistance.lineEdit.text()))
        dic = IP.folowMaxValue(image2Folow, seed1, seed2, maxDistance)

        diameter = dic['Diameter']
        distance = dic['Distance']
        Image = dic['Image']
        np.savetxt('./Data/Line.txt', [diameter, distance])

        self.plt.show()
        self.plt.clear()
        self.plt.showGrid(x=True, y=True)
        self.plt.setLabel('left', 'Diameter')
        self.plt.setLabel('bottom', 'Distance')
        self.curve = self.plt.plot(distance, diameter)
        self.plt.repaint()
        self.addImage('FollowedLine', Image)

    def stack_Function(self):

        z11 = int(str(self.index_z1_i1.lineEdit.text()))
        z21 = int(str(self.index_z2_i1.lineEdit.text()))
        z12 = int(str(self.index_z1_i2.lineEdit.text()))
        z22 = int(str(self.index_z2_i2.lineEdit.text()))

        im1 = self.Data_list[self.Im1.currentIndex()]
        im2 = self.Data_list[self.Im2.currentIndex()]

        if (im1.shape[1] == im2.shape[1]) and (im1.shape[2] == im2.shape[2]):
            if (z21 > z11) and (z22 > z12):

                ImageOut = np.zeros(((z21 - z11) + (z22 - z12), im1.shape[1], im1.shape[2]))
                ImageOut[0:(z21 - z11), :, :] = im1[z11:z21, :, :]
                ImageOut[(z21 - z11):(z21 - z11) + (z22 - z12), :, :] = im2[z12:z22, :, :]
                self.addImage("Stack_ ", ImageOut)
            else:
                msgBox = qt.QMessageBox(self)
                msgBox.setText('Unconsistent coordinate ')
                msgBox.exec_()

        else:
            msgBox = qt.QMessageBox(self)
            msgBox.setText('Images Size don t corespond ')
            msgBox.exec_()

        self._stack_GUI()

    def destack_Function(self):

        h = int(self.index_h.lineEdit.text())
        zmin = int(self.index_zmin.lineEdit.text())
        zmax = int(self.index_zmax.lineEdit.text())

        im1 = self.Data_list[self.imageSelection.currentRow()]
        stageNumber = int(im1.shape[0] / h)

        for z in range(int(h / 2.0), stageNumber * h, h):
            if ((z - zmin) > 0) and ((zmax + z) < im1.shape[0]):
                hZ = int(h / 2.0)
                stack = im1[z - (hZ - zmin):z + (zmax - hZ), :]
                self.addImage(str(z), stack, '', px_size)

        self._destack_GUI()

    def hessian_GUI(self):
        self.hide_all_button()

        image = self.Data_list[self.imageSelection.currentRow()]
        H = IP.hessian_matrix(image)
        im1, im2, im3 = IP.computeEigenvalueHessianMatrix(H)


        self.addImage("Hessian1_ ",im1)
        self.addImage("Hessian2_ ",im2)
        self.addImage("Hessian3_ ",im3)

    def frangi_GUI(self):
        self.hide_all_button()
        self.MIButton = qt.QPushButton("Frangi !")
        self.scaleRangeMin = LabelEditAndButton(True, "scale Range Min: ", True, str(0), False)
        self.scaleRangeMax = LabelEditAndButton(True, "scale Range Min: ", True, str(10), False)
        self.scale_step = LabelEditAndButton(True, "scale_step: ", True, str(2), False)
        self.alpha = LabelEditAndButton(True, "Alpha: ", True, str(0.5), False)
        self.beta = LabelEditAndButton(True, "Beta: ", True, str(0.5), False)
        self.frangi_c = LabelEditAndButton(True, "Frangi Coefficient: ", True, str(10), False)
        self.BlackVessels = qt.QCheckBox("Black Vessels")

        self.MIButton.setMaximumWidth(250)
        self.scaleRangeMin.setMaximumWidth(250)
        self.scaleRangeMax.setMaximumWidth(250)
        self.scale_step.setMaximumWidth(250)
        self.alpha.setMaximumWidth(250)
        self.beta.setMaximumWidth(250)
        self.frangi_c.setMaximumWidth(250)
        self.BlackVessels.setMaximumWidth(250)

        self.MIButton.clicked.connect(self.frangi)

        self.buttonLayout.addWidget(self.scaleRangeMin)
        self.buttonLayout.addWidget( self.scaleRangeMax)
        self.buttonLayout.addWidget(self.scale_step)
        self.buttonLayout.addWidget(self.alpha)
        self.buttonLayout.addWidget(self.beta)
        self.buttonLayout.addWidget(self.frangi_c)
        self.buttonLayout.addWidget( self.BlackVessels)
        self.buttonLayout.addWidget(self.MIButton)


    def frangi(self):

        image = self.Data_list[self.imageSelection.currentRow()]

        scale_range_min = int(self.scaleRangeMin.lineEdit.text())
        scale_range_max = int(self.scaleRangeMax.lineEdit.text())
        scale_step =  int(self.scale_step.lineEdit.text())
        alpha = float(self.beta.lineEdit.text())
        beta =  float(self.beta.lineEdit.text())
        frangi_c = float(self.frangi_c.lineEdit.text())

        if self.BlackVessels.checkState() == 2 :
            Vessels = True
        else:
            Vessels = False

        Image = IP.frangi(image,  (scale_range_min, scale_range_max),scale_step, alpha,beta, frangi_c,Vessels)

        self.addImage("Frangi", Image )


    def _analysisROIGUI(self):

        self.hide_all_button()
        if any(self.ItemsLists[self.imageSelection.currentRow()]["Zones"]["Direction0"]):

            Zone = self.ItemsLists[self.imageSelection.currentRow()]["Zones"]["Direction0"][0]

            ROI = [Zone[2], Zone[0], Zone[3], Zone[1], Zone[4]]

            roi_2_analyse = self.Data_list[self.imageSelection.currentRow()][ROI[0], ROI[1]:ROI[2], ROI[3]:ROI[4]]

            nb_px = (ROI[2] - ROI[1]) * (ROI[4] - ROI[3])

            mean_v = np.mean(roi_2_analyse)
            std_v = np.var(roi_2_analyse)

        elif any(self.ItemsLists[self.imageSelection.currentRow()]["Circles"]["Direction0"]):

            Circle = self.ItemsLists[self.imageSelection.currentRow()]["Circles"]["Direction0"][0]
            nb_px = np.pi * (Circle[3] / 2.0) ** 2.0
            imageToAnalyse = self.Data_list[self.imageSelection.currentRow()][Circle[2], :, :]

            a = Circle[1]
            b = Circle[0]

            nx, ny = imageToAnalyse.shape
            y, x = np.ogrid[-a:nx - a, -b:ny - b]
            mask = x * x + y * y <= ((Circle[3] / 2.0) ** 2.0)
            mean_v = np.mean(imageToAnalyse[mask])
            std_v = np.sqrt(np.var(imageToAnalyse[mask]))
            nb_px = np.count_nonzero(mask)

        text = "Distribution of the selected ROI\n"
        text += 'Number of Pixels : ' + str(
            nb_px) + '\n' + 'Mean Value : ' + '%.10f' % mean_v + '\n' + 'Standart Deviation : ' + '%.10f' % std_v + '\n'
        self.hide_all_button()
        self.txtDisplay = qt.QTextEdit("Results")
        self.txtDisplay.setText(text)
        self.txtDisplay.setMaximumWidth(250)
        self.buttonLayout.addWidget(self.txtDisplay)

    def _MIGUI(self):
        self.hide_all_button()
        self.MIButton = qt.QPushButton("Compute Mutual Information")
        self.startFrame = LabelEditAndButton(True, "First Frame: ", True, str(0), False)
        self.endFrame = LabelEditAndButton(True, "Last Frame: ", True, str(100), False)

        self.MIButton.setMaximumWidth(250)
        self.startFrame.setMaximumWidth(250)
        self.endFrame.setMaximumWidth(250)

        self.MIButton.clicked.connect(self.MI_Function)

        self.buttonLayout.addWidget(self.startFrame)
        self.buttonLayout.addWidget(self.endFrame)
        self.buttonLayout.addWidget(self.MIButton)

    def MI_Function(self):
        Image = self.Data_list[self.imageSelection.currentRow()]

        if len(self.ItemsLists[self.imageSelection.currentRow()]["Zones"]["Direction0"]) != 0:
            FlagROI = True
            Zone = self.ItemsLists[self.imageSelection.currentRow()]["Zones"]["Direction0"][0]
            ROI = [Zone[0], Zone[3], Zone[1], Zone[4]]

        else:
            FlagROI = False

        startFrame = int(self.startFrame.lineEdit.text())
        endFrame = int(self.endFrame.lineEdit.text())

        arrayMI = IP.mutualInformation(Image, [startFrame, endFrame], FlagROI, ROI)

        array_time = np.arange(startFrame, endFrame - 1)

        self.plt.show()
        self.plt.clear()
        self.plt.showGrid(x=True, y=True)
        self.plt.setLabel('left', 'Mutual Information')
        self.plt.setLabel('bottom', 'Image Numbers')
        self.curve = self.plt.plot(array_time, arrayMI)
        self.plt.repaint()

    def correct_ff_GUI(self):

        self.hide_all_button()
        self.corButton = qt.QPushButton("Correct")
        self.label1 = qt.QLabel("Projection Images")
        self.ImP = qt.QComboBox()
        self.label2 = qt.QLabel("Flat Field Image")
        self.ImFF = qt.QComboBox()

        self.corButton.setMaximumWidth(250)
        self.label1.setMaximumWidth(250)
        self.ImP.setMaximumWidth(250)
        self.label2.setMaximumWidth(250)
        self.ImFF.setMaximumWidth(250)

        self.ImP.addItems(self.Name_list)
        self.ImFF.addItems(self.Name_list)

        self.corButton.connect.clicked(self.correct_ff_function)

        self.buttonLayout.addWidget(self.label1)
        self.buttonLayout.addWidget(self.ImP)
        self.buttonLayout.addWidget(self.label2)
        self.buttonLayout.addWidget(self.ImFF)
        self.buttonLayout.addWidget(self.corButton)

    def correct_ff_function(self):
        Prj = self.Data_list[self.ImP.currentIndex()]
        FF = self.Data_list[self.ImFF.currentIndex()]
        FFMed = np.copy(np.median(FF, axis=0) + 1.0)

        Norm = np.zeros((Prj.shape[0], Prj.shape[1], Prj.shape[2]))
        for z in range(0, Prj.shape[0]):
            Norm[z, :, :] = Prj[z, :, :] / (FFMed)

        self.addImage('Norm_' + self.Name_list[self.ImP.currentIndex()], Norm)

    def _FFTCORGUI(self):

        self.hide_all_button()

        self.corButton = qt.QPushButton("Correlate")
        self.label1 = qt.QLabel("Fixed Image")
        self.ImF = qt.QComboBox()
        self.sizeMaskF = LabelEditAndButton(True, "Size Fixed Image Mask: ", True, str(100), False)

        self.label2 = qt.QLabel("Moving Image")
        self.ImM = qt.QComboBox()
        self.sizeMaskM = LabelEditAndButton(True, "Size Fixed Image Mask: ", True, str(100), False)

        self.overlap = LabelEditAndButton(True, "Overlap Ratio (0-1): ", True, str(0.0), False)

        self.corButton.setMaximumWidth(250)
        self.label1.setMaximumWidth(250)
        self.ImF.setMaximumWidth(250)
        self.sizeMaskF.setMaximumWidth(250)
        self.label2.setMaximumWidth(250)
        self.ImM.setMaximumWidth(250)
        self.sizeMaskM.setMaximumWidth(250)

        self.overlap.setMaximumWidth(250)

        self.ImF.addItems(self.Name_list)
        self.ImM.addItems(self.Name_list)

        self.corButton.clicked.connect(self.CORFFT_Function)

        self.buttonLayout.addWidget(self.label1)
        self.buttonLayout.addWidget(self.ImF)
        self.buttonLayout.addWidget(self.sizeMaskF)
        self.buttonLayout.addWidget(self.label2)
        self.buttonLayout.addWidget(self.ImM)
        self.buttonLayout.addWidget(self.sizeMaskM)
        self.buttonLayout.addWidget(self.overlap)

        self.buttonLayout.addWidget(self.corButton)

    def CORFFT_Function(self):

        sizeF = int(str(self.sizeMaskF.lineEdit.text()))
        sizeM = int(str(self.sizeMaskM.lineEdit.text()))

        overlap = float(str(self.overlap.lineEdit.text()))
        ImF = self.Data_list[self.ImF.currentIndex()]
        ImM = self.Data_list[self.ImM.currentIndex()]

        ImageOut = IP.FFTCOR(ImF, ImM, sizeF, sizeM, overlap)

        self.addImage("Corr ", ImageOut)
        self._FFTCORGUI()

    def _RatioMaskGUI(self):

        self.hide_all_button()

        self.ratioButton = qt.QPushButton("Start")

        self.WindowSize = LabelEditAndButton(True, "Window Size", True, "4", False)
        self.Overlap = LabelEditAndButton(True, "Overlap", True, "0.5", False)
        self.ratioButton.clicked.connect(self.study_maskRatio)

        self.ratioButton.setMaximumWidth(250)
        self.WindowSize.setMaximumWidth(250)
        self.Overlap.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.WindowSize)
        self.buttonLayout.addWidget(self.Overlap)
        self.buttonLayout.addWidget(self.ratioButton)

    def _RatioVolMaskGUI(self):

        self.hide_all_button()

        self.ratioButton = qt.QPushButton("Start")

        self.WindowSize = LabelEditAndButton(True, "Window Size", True, "4", False)
        self.Overlap = LabelEditAndButton(True, "Overlap", True, "0.5", False)
        self.ratioButton.clicked.connect(self.study_maskVolRatio)

        self.ratioButton.setMaximumWidth(250)
        self.WindowSize.setMaximumWidth(250)
        self.Overlap.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.WindowSize)
        self.buttonLayout.addWidget(self.Overlap)
        self.buttonLayout.addWidget(self.ratioButton)

    def study_maskRatio(self):
        Image = self.Data_list[self.imageSelection.currentRow()]

        ws = int(self.WindowSize.lineEdit.text())
        ov = float(self.Overlap.lineEdit.text())

        DensityMap = IP.ComputeMaskLocalDensity(np.copy(Image), ws, ov)

        self.addImage("DensityNumberMap", DensityMap)

    def study_maskVolRatio(self):
        Image = self.Data_list[self.imageSelection.currentRow()]

        ws = int(self.WindowSize.lineEdit.text())
        ov = float(self.Overlap.lineEdit.text())

        DensityMap = IP.ComputeMaskLocalVolumeDensity(np.copy(Image), ws, ov)

        self.addImage("DensityVolumeMap", DensityMap)

    """ Morphologie """

    def _morpho_GUI(self):

        self.hide_all_button()

        self.morphoButton = qt.QPushButton("Compute")
        self.Image = qt.QComboBox()
        self.checkErosion = qt.QCheckBox("Erosion")
        self.checkDilatation = qt.QCheckBox("Dilatation")
        self.checkOpening = qt.QCheckBox("Opening")
        self.checkClosing = qt.QCheckBox("Closing")
        self.checkFilling = qt.QCheckBox("Filling 2D")
        self.checkDistance = qt.QCheckBox("Distance")
        self.checkSkeletton = qt.QCheckBox("Skeletton")
        self.kernelRadius = LabelEditAndButton(True, "Kernel radius", True, "1", False)

        self.Image.addItems(self.Name_list)
        self.kernelRadius.setFixedSize(250, 40)

        self.morphoButton.clicked.connect(self.morphoFunction)
        self.checkErosion.stateChanged.connect(self.Ero)
        self.checkDilatation.stateChanged.connect(self.Dil)
        self.checkOpening.stateChanged.connect(self.Open)
        self.checkClosing.stateChanged.connect(self.Closing)
        self.checkFilling.stateChanged.connect(self.Filling)
        self.checkDistance.stateChanged.connect(self.Distance)
        self.checkSkeletton.stateChanged.connect(self.Skeletton)

        self.buttonLayout.addWidget(self.Image)
        self.buttonLayout.addWidget(self.checkErosion)
        self.buttonLayout.addWidget(self.checkDilatation)
        self.buttonLayout.addWidget(self.checkOpening)
        self.buttonLayout.addWidget(self.checkClosing)
        self.buttonLayout.addWidget(self.checkFilling)
        self.buttonLayout.addWidget(self.checkDistance)
        self.buttonLayout.addWidget(self.checkSkeletton)
        self.buttonLayout.addWidget(self.kernelRadius)
        self.buttonLayout.addWidget(self.morphoButton)

    def morphoFunction(self):

        npArray = self.Data_list[self.Image.currentIndex()]

        kernel = int(self.kernelRadius.lineEdit.text())

        if self.checkDilatation.checkState() == 2:
            name = 'Dil' + str(kernel)
            ImageOut = IP.morpho('Dilate', npArray, kernel)

        if self.checkErosion.checkState() == 2:
            name = 'Ero' + str(kernel)
            ImageOut = IP.morpho('Erode', npArray, kernel)

        if self.checkOpening.checkState() == 2:
            name = 'Open' + str(kernel)
            ImageOut = IP.morpho('Open', npArray, kernel)

        if self.checkClosing.checkState() == 2:
            name = 'Clos' + str(kernel)
            ImageOut = IP.morpho('Close', npArray, kernel)

        if self.checkFilling.checkState() == 2:
            name = 'Fill' + str(kernel)
            ImageOut = IP.morpho('Fill', npArray, kernel)

        if self.checkDistance.checkState() == 2:
            name = 'Distance' + str(kernel)
            ImageOut = IP.Distance(npArray)

        if self.checkSkeletton.checkState() == 2:
            name = 'Skeletton' + str(kernel)
            ImageOut = IP.morpho('Skeletton', npArray, kernel)

        self.addImage(name, ImageOut)

        self._morpho_GUI()

    """ Registration """

    def _registration_GUI(self):
        self.optionRegisterW = RegisteringOption(self.Name_list)
        self.optionRegisterW.startRegistering.clicked.connect(self.lauchRegistering)
        self.optionRegisterW.show()

    def lauchRegistering(self, dicPar=-1):
        if dicPar == False:
            dicPar = self.optionRegisterW.dicPar

        self.plt.show()
        self.plt.clear()
        pen_act = pg.mkPen((0, 200, 0, 255), width=1)
        pen_des = pg.mkPen((200, 0, 0, 255), width=1)
        pen_bot = pg.mkPen((200, 0, 200, 255), width=1)

        self.p1 = self.plt.addPlot(row=0, col=0, title="Metric Time Lapse ",
                                   labels={'left': "Log Metric Value ", 'bottom': "Iteration Number"})
        self.p2 = self.plt.addPlot(row=1, col=0, title="Convergence Time Lapse",
                                   labels={'left': "Log Convergence", 'bottom': "Iteration Number"})
        self.p3 = self.plt.addPlot(row=2, col=0, title="Learning Rate Time Lapse",
                                   labels={'left': "Log Learning Rate", 'bottom': "Iteration Number"})

        self.curve1 = self.p1.plot([], pen=pen_act)
        self.curve2 = self.p2.plot([], pen=pen_des)
        self.curve3 = self.p3.plot([], pen=pen_bot)

        # self.curve1 = self.p1.plot()

        #        self.plt.clear()
        #        self.plt.showGrid(x=True, y=True)
        #        self.plt.setLabel('left', 'Log(Metric)')
        #        self.plt.setLabel('bottom', 'Iterations')
        #        self.curve = self.plt.plot()
        #        self.curve.setData( [])
        #        self.plt.repaint()

        self.counterIter = 0

        self.R = IPR.Registering(dicPar, self.Data_list, [1, 1, 1])
        self.regThread = rT.registerThread(self.R, self)

        self.regThread.RegDone.connect(self.endRegister)

        self.R.reg_method.AddCommand(sitk.sitkStartEvent, self.start_plot)
        self.R.reg_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, self.update_multires_iterations)
        self.R.reg_method.AddCommand(sitk.sitkIterationEvent, lambda: self.plot_values())

        self.iterNumber = 0
        self.start_time = time.time()
        self.regThread.start()

    """Utility"""

    def start_plot(self):
        self.metric_values = []
        self.converg_values = []
        self.lr_values = []
        self.multires_iterations = []

    def endRegister(self):

        print("--- %s seconds ---" % (time.time() - self.start_time))
        newImage = self.R.returnNumpyImage()
        self.addImage("IMovingBefore", self.R.ImageMovingBef)
        self.addImage("IFixedBefore", self.R.ImageFixedBef)
        ImageMoving = np.copy(self.R.ImageMovingBef)
        ImageFixe = np.copy(self.R.ImageFixedBef)
        self.addImage("RegiterIm", newImage)
        self.addImage("Jacobian", self.R.vectorFieldJacobian)

        MapChess = np.copy(IP.giveChessImage(newImage, ImageFixe, 50))
        self.addImage("Chess", MapChess)
        vectorField = np.copy(self.R.transformOut)

        normVectorField = np.copy(np.sqrt(
            np.square(vectorField[:, :, :, 0] ) + np.square(vectorField[:, :, :, 1] ) + np.square(
                vectorField[:, :, :, 2] )))

        self.addImage("Transform", normVectorField)

        size_x = normVectorField.shape[0]
        size_y = normVectorField.shape[1]
        size_z = normVectorField.shape[2]

        np.save('./VectorField' + self.Name_list[0], vectorField)

        step = 5

        for x in range(0, size_x - 1, step):
            for y in range(0, size_y - 1, step):
                for z in range(0, size_z - 1, step):

                    x0 = x
                    y0 = y
                    z0 = z
                    x1 = vectorField[x, y, z, 0] + x
                    y1 = vectorField[x, y, z, 1] + y
                    z1 = vectorField[x, y, z, 2] + z

                    if (int(x1) < ImageFixe.shape[0]) and (int(y1) < ImageFixe.shape[1]) and (
                            int(z1) < ImageFixe.shape[2]):
                        if (ImageFixe[int(x1), int(y1), int(z1)] != 0):
                            if (x < ImageMoving.shape[0]) and (y < ImageMoving.shape[1]) and (z < ImageMoving.shape[2]):
                                if (ImageMoving[x, y, z] != 0):
                                    self.ItemsLists[self.imageSelection.currentRow() - 4]['Arrows'][
                                        'Direction0'].append([z0, y0, x0, z1, y1, x1])
                                    self.ItemsLists[self.imageSelection.currentRow() - 4]['Arrows'][
                                        'Direction1'].append([x0, y0, z0, x1, y1, z1])
                                    self.ItemsLists[self.imageSelection.currentRow() - 4]['Arrows'][
                                        'Direction2'].append([x0, y0, z0, x1, y1, z1])
                                    self.ItemsLists[self.imageSelection.currentRow() - 3]['Arrows'][
                                        'Direction0'].append([z0, y0, x0, z1, y1, x1])
                                    self.ItemsLists[self.imageSelection.currentRow() - 3]['Arrows'][
                                        'Direction1'].append([x0, y0, z0, x1, y1, z1])
                                    self.ItemsLists[self.imageSelection.currentRow() - 3]['Arrows'][
                                        'Direction2'].append([x0, y0, z0, x1, y1, z1])
                                    self.ItemsLists[self.imageSelection.currentRow()]['Arrows']['Direction0'].append(
                                        [z0, y0, x0, z1, y1, x1])
                                    self.ItemsLists[self.imageSelection.currentRow()]['Arrows']['Direction1'].append(
                                        [x0, y0, z0, x1, y1, z1])
                                    self.ItemsLists[self.imageSelection.currentRow()]['Arrows']['Direction2'].append(
                                        [x0, y0, z0, x1, y1, z1])

    #    def change_GUI_Register(self):
    #        if self.iterNumber == 0:
    #            image = self.R.ApplyMovingTransform([150,150,150])
    #            print image[0].shape[0]
    #            print image[0].shape[1]
    #            self.track_image = np.zeros((500,image[0].shape[0],image[0].shape[1]))
    #            self.track_image[0,:,:] = image[0]
    #            self.addImage('Test',self.track_image)
    #        elif self.iterNumber < 500:
    #            self.track_image[self.iterNumber,:,:]= self.R.ApplyMovingTransform([150,150,150])[0]
    #            self.image3DWidget.axialWidget.sliceSlider.slider.setTickPosition(self.iterNumber)
    #            self.image3DWidget.axialWidget._changeSlice()
    #
    #        self.curve1.setData(self.metric_values, pen= pg.mkPen((0,200,0,255), width=1))
    #        self.curve2.setData(self.converg_values, pen=pg.mkPen((200,0,0,255), width=1))
    #        self.curve3.setData(self.lr_values, pen=pg.mkPen((200,0,200,255), width=1))
    #
    #        print 'IN'

    def RegisterIter(self):
        print('InIter')

    def plot_values(self):

        # ThreadDisplay = rT.registerDisplayThread(self)
        # ThreadDisplay.start()

        MetricValue = self.R.reg_method.GetMetricValue()
        # CurrentLevel = self.R.reg_method.GetCurrentLevel()
        ConvergenceValue = self.R.reg_method.GetOptimizerConvergenceValue()
        LearningRate = self.R.reg_method.GetOptimizerLearningRate()

        self.metric_values.append(np.log(MetricValue))

        if ConvergenceValue < 100000.0:
            self.converg_values.append(np.log(ConvergenceValue))
        else:
            self.converg_values.append(np.log(100000.0))

        if LearningRate < 100000.0:
            self.lr_values.append(np.log(LearningRate))
        else:
            self.lr_values.append(np.log(100000.0))

        self.curve1.setData(self.metric_values, pen=pg.mkPen((0, 200, 0, 255), width=1))
        self.curve2.setData(self.converg_values, pen=pg.mkPen((200, 0, 0, 255), width=1))
        self.curve3.setData(self.lr_values, pen=pg.mkPen((200, 0, 200, 255), width=1))
        # OptimizerPosition = self.R.reg_method.GetOptimizerPosition()
        # OptimizerScales = self.R.reg_method.GetOptimizerScales()
        # StopConditionDescription = self.R.reg_method.GetOptimizerStopConditionDescription()

        # self.emit(qt.SIGNAL('self.change_GUI_Register()'))
        # Transform = self.R.initial_transform
        # self.R.ApplyTransform(central_indexes)

        # print '------------------------------------'
        # print 'CurrentLevel :',CurrentLevel

        # print 'LearningRate :', LearningRate
        # 'OptimizerPosition :',OptimizerPosition
        # print 'OptimizerScales :', OptimizerScales
        # print 'StopConditionDescription : ',StopConditionDescription
        # print 'Transfrom :', Transform

        # fill = pg.FillBetweenItem(self.curve1, self.curve2, brush=(100,0,100,100))
        # self.p1.addItem(fill)
        # self.curve.setData(np.log(self.metric_values))
        self.iterNumber += 1

    def update_multires_iterations(self):
        self.counterIter += 1

    def ExpiCheck(self):
        if self.Expi.checkState() == 2:
            self.Inspi.setCheckState(0)

    def InspiCheck(self):
        if self.Inspi.checkState() == 2:
            self.Expi.setCheckState(0)

    """ GUI Study """

    def changeROIMeasure(self):
        gap = 3.0

        self.z_value1_t = int(self.pos1.lineEdit.text()) + (gap / 2.0)
        self.z_value1_b = int(self.pos1.lineEdit.text()) - (gap / 2.0)
        self.z_value2_t = int(self.pos2.lineEdit.text()) + (gap / 2.0)
        self.z_value2_b = int(self.pos2.lineEdit.text()) - (gap / 2.0)
        self.z_value3_t = int(self.pos3.lineEdit.text()) + (gap / 2.0)
        self.z_value3_b = int(self.pos3.lineEdit.text()) - (gap / 2.0)

        x_value_2 = self.image3DWidget.sagittalWidget.sliceSlider.value()
        x_value_1 = self.image3DWidget.coronalWidget.sliceSlider.value()

        y_value_1_b = 0
        y_value_1_t = self.Data_list[self.imageSelection.currentRow()].shape[1]
        y_value_2_b = 0
        y_value_2_t = self.Data_list[self.imageSelection.currentRow()].shape[2]

        self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction1'] = []
        self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction2'] = []

        for x in range(0, self.Data_list[self.imageSelection.currentRow()].shape[2]):
            self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction1'].append(
                [y_value_1_b, x, self.z_value1_b, y_value_1_t, x, self.z_value1_t])
            self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction1'].append(
                [y_value_1_b, x, self.z_value2_b, y_value_1_t, x, self.z_value2_t])
            self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction1'].append(
                [y_value_1_b, x, self.z_value3_b, y_value_1_t, x, self.z_value3_t])

        for x in range(0, self.Data_list[self.imageSelection.currentRow()].shape[1]):
            self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction2'].append(
                [x, y_value_2_b, self.z_value1_b, x, y_value_2_t, self.z_value1_t])
            self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction2'].append(
                [x, y_value_2_b, self.z_value2_b, x, y_value_2_t, self.z_value2_t])
            self.ItemsLists[self.imageSelection.currentRow()]['Zones']['Direction2'].append(
                [x, y_value_2_b, self.z_value3_b, x, y_value_2_t, self.z_value3_t])

        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])


    def endSegmentationReg(self):

        self.hide_all_button()

        self.save = qt.QPushButton()
        self.save.setIcon(qt.QIcon(qt.QPixmap('./Icones/save.png')))
        self.save.setMaximumWidth(250)
        self.save.clicked.connect(self._buttonSaveSeg)

        self.restartProcess = qt.QPushButton('Restart Segmentation')
        self.restartProcess.clicked.connect(self._buttonRestartSeg)

        self.stepSize = LabelEditAndButton(True, "Step Size", True, "10.0", False)
        self.IterationNumber = LabelEditAndButton(True, "Iteration Number", True, "100", False)
        self.ResolutionReg = LabelEditAndButton(True, "Resolution Reg", True, "4", False)

        self.stepSize.setMaximumWidth(250)
        self.IterationNumber.setMaximumWidth(250)
        self.ResolutionReg.setMaximumWidth(250)

        self.startProcess = qt.QPushButton('StartRegistration')
        self.startProcess.clicked.connect(self._buttonStartRegistration)

        self.stepSize.setMaximumWidth(250)

        self.buttonLayout.addWidget(self.save)
        self.buttonLayout.addWidget(self.restartProcess)
        self.buttonLayout.addWidget(self.stepSize)
        self.buttonLayout.addWidget(self.IterationNumber)
        self.buttonLayout.addWidget(self.ResolutionReg)
        self.buttonLayout.addWidget(self.startProcess)

    def _buttonStartRegistration(self):

        IndexIFIT = self.Seg3Index + 1
        IndexMI = self.Seg2Index
        IndexIMIT = self.Seg4Index + 1
        IndexFI = self.Seg1Index

        TracheaVolI = np.copy(self.Data_list[IndexIFIT - 1])

        TracheaVolI[TracheaVolI == 1] = -1
        TracheaVolI += 1

        ImageExpiR = np.copy(self.Data_list[IndexMI])
        ImageExpiR *= TracheaVolI

        mask = np.ones((ImageExpiR.shape[0], ImageExpiR.shape[1], ImageExpiR.shape[1]))
        mask[ImageExpiR < -1100.0] = 0
        mask[ImageExpiR > -500] = 0
        self.volumeExpi = np.sum(np.sum(mask))

        ImageInspiR = self.Data_list[IndexFI]
        ImageInspiR *= TracheaVolI

        mask = np.ones((ImageInspiR.shape[0], ImageInspiR.shape[1], ImageInspiR.shape[1]))
        mask[ImageInspiR < -1100.0] = 0
        mask[ImageInspiR > -500] = 0
        self.volumeInspi = np.sum(np.sum(mask))

        stepSize = float(str(self.stepSize.lineEdit.text()))
        Iter = int(str(self.IterationNumber.lineEdit.text()))
        Res = int(str(self.ResolutionReg.lineEdit.text()))

        if Res == 16:
            Scaling = [16, 16, 16, 16, 0, 0, 0, 0]
        elif Res == 8:
            Scaling = [16, 16, 8, 8, 0, 0, 0, 0]
        elif Res == 4:
            Scaling = [16, 8, 8, 4, 0, 0, 0, 0]
        elif Res == 2:
            Scaling = [16, 8, 4, 2, 0, 0, 0, 0]

        nNodesX = 9
        nNodesY = 9
        nNodesZ = 10

        shapeX = self.Data_list[self.Seg1Index].shape[1]
        shapeY = self.Data_list[self.Seg1Index].shape[2]
        shapeZ = self.Data_list[self.Seg1Index].shape[0]

        GridX = int(shapeX / (nNodesX - 2))
        GridY = int(shapeY / (nNodesY - 2))
        GridZ = int(shapeZ / (nNodesZ - 2))

        dicPar = {'Inputs': {'MIM': 0, 'IFIT': IndexIFIT, 'MI': IndexMI, 'IMIT': IndexIMIT, 'InitT': u'Moments',
                            'FI': IndexFI, 'FIM': 0}, \
                  'Optimizer': {'Par': [stepSize, Iter, 1e-07, 20, 0.01, stepSize, 0.05, 10.0, 2.0, 10.0],
                                'ScalePar': [5.0, 0.01], \
                                'Method': 'Conjugate Gradient Line Search', 'MethodScaling': 'Index Shift'}, \
                  'Outputs': {'Save': [], 'Display': [], 'Sampling': 2.0}, \
                  'Metric': {'GradM': 1, 'Par': [], 'GradF': 1, 'Method': 'Means Squares', \
                             'Sampling': {'Percentage': 0.5, 'Method': u'Random'}}, 'Scaling': Scaling, \
                  'Interpolator': 'BSpline', 'Grid': [GridX, GridY, GridZ]}

        self.lauchRegistering(dicPar)


    def addImage(self, name, Image, tooltip="", pxSize=[1, 1, 1]):

        item = qt.QListWidgetItem(name)
        item.setToolTip(tooltip)
        item.setFlags(item.flags() | qt.Qt.ItemIsEditable)
        self.imageSelection.itemChanged.connect(self.nameImageChange)
        self.Name_list.append(name)
        item.setTextAlignment(qt.Qt.AlignHCenter)
        self.imageSelection.addItem(item)

        self.Data_list.append(Image)
        self.imageSelection.setCurrentRow(self.imageSelection.count() - 1)
        self.ItemsLists.append(copy.deepcopy(self.ItemsInit))
        self.Overlays.append(copy.deepcopy(self.OverlayPar))

        self.image3DWidget.updateItems(self.ItemsLists[self.imageSelection.currentRow()])
        self.image3DWidget.updateOverlays(self.Overlays[self.imageSelection.currentRow()])

    def hide_all_button(self):
        self.plt.hide()

        for i in range(self.buttonLayout.count()):
            if i != 0:
                self.buttonLayout.itemAt(i).widget().hide()
