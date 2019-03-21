# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:26:20 2016

@author: broche
"""
import sys
import PyMcaQt as qt
from Frame import Frame


class MyMainWindow(qt.QMainWindow) :
    def __init__(self, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.mainWidget = Frame(self)
        self.setCentralWidget(self.mainWidget)
        self.setWindowTitle('Pini3_Alpha')
        self._buildMenu()

    def _buildMenu(self) :
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        loadAction = qt.QAction('&Load Slices Image', self)
        loadAction.setShortcut('Ctrl+O')
        loadAction.setStatusTip('Open 3D Images')
        loadAction.triggered.connect(self.mainWidget._load)
        fileMenu.addAction(loadAction)

        loadActionF = qt.QAction('&Load Folders Image', self)
        loadActionF.setStatusTip('Open 3D Images')
        loadActionF.triggered.connect( self.mainWidget._loadFolders)
        fileMenu.addAction(loadActionF)


        loadDicomAction = qt.QAction('&Load Dicom/nftii/Mat File', self)
        loadDicomAction.setStatusTip('Load Dicom/Mat File')
        loadDicomAction.triggered.connect(self.mainWidget._loadDicom)
        fileMenu.addAction(loadDicomAction)

        loadSTLAction = qt.QAction('&Load STL file as Image', self)
        loadSTLAction.setStatusTip('Load STL File')
        loadSTLAction.triggered.connect(self.mainWidget._loadSTL)
        fileMenu.addAction(loadSTLAction)

        saveAction = qt.QAction( '&Save', self)
        saveAction.setShortcut('Ctrl+S')
        saveAction.setStatusTip('Save ')
        saveAction.triggered.connect(self.mainWidget._save)
        fileMenu.addAction(saveAction)

        exitAction = qt.QAction(qt.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(qt.qApp.quit)
        fileMenu.addAction(exitAction)

        ProcessMenu = menubar.addMenu('&Image Processing')

        FormatingMenu = ProcessMenu.addMenu('&Image Format')
        cropAction = qt.QAction('&Crop', self)
        cropAction.setStatusTip('Crop Volume')
        cropAction.triggered.connect( self.mainWidget._crop )
        FormatingMenu.addAction(cropAction)

        fillAction = qt.QAction('&Fill', self)
        fillAction.setStatusTip('Fill Volume')
        fillAction.triggered.connect( self.mainWidget._fill )
        FormatingMenu.addAction(fillAction)

        resampleAction = qt.QAction('&Resample', self)
        resampleAction.triggered.connect( self.mainWidget._resample_GUI )
        FormatingMenu.addAction(resampleAction)

        contrastAction = qt.QAction('&Contrast', self)
        contrastAction.triggered.connect( self.mainWidget._contrast_GUI )
        FormatingMenu.addAction(contrastAction)

        rotationAction = qt.QAction('&Rotation', self)
        rotationAction.triggered.connect(self.mainWidget._rotation_GUI )
        FormatingMenu.addAction(rotationAction)

        equalizeAction = qt.QAction('&Equalize', self)
        equalizeAction.setStatusTip('Equalize Slices using 1 or 2 uniform zone')
        equalizeAction.triggered.connect( self.mainWidget._equalize)
        FormatingMenu.addAction(equalizeAction)

        equalizeActionPath = qt.QAction('&Equalize along Path', self)
        equalizeActionPath.setStatusTip('Equalize along path')
        equalizeActionPath.triggered.connect( self.mainWidget._equalizePath)
        FormatingMenu.addAction(equalizeActionPath)

        equalizeHistoAction = qt.QAction('&Equalize Histogram', self)
        equalizeHistoAction.setStatusTip('Put the 2 Main Modes of the histogram to a specific Value')
        equalizeHistoAction.triggered.connect(self.mainWidget._equalizeHisto)
        FormatingMenu.addAction(equalizeHistoAction)

        mozaicAction = qt.QAction('&Mozaic', self)
        mozaicAction.setStatusTip('Make a Mozaic Stack')
        mozaicAction.triggered.connect(self.mainWidget._mozaicGUI)
        FormatingMenu.addAction(mozaicAction)

        stackAction = qt.QAction('&Stack Images', self)
        stackAction.setStatusTip('Stack Images Together')
        stackAction.triggered.connect(self.mainWidget._stack_GUI)
        FormatingMenu.addAction(stackAction)

        destackAction = qt.QAction('&Destack Images', self)
        destackAction.setStatusTip('Destack Images Together')
        destackAction.triggered.connect(self.mainWidget._destack_GUI)
        FormatingMenu.addAction(destackAction)

        AnalyseMenu = ProcessMenu.addMenu('&Analyse')
        analyseAction = qt.QAction('&Analyse Distribution In ROI', self)
        analyseAction.triggered.connect( self.mainWidget._analysisROIGUI)
        AnalyseMenu.addAction(analyseAction)

        analyseLineAction = qt.QAction('&Read Pixels Line', self)
        analyseLineAction.triggered.connect( self.mainWidget._followLineGUI)
        AnalyseMenu.addAction( analyseLineAction)

        analyseMaxLineAction = qt.QAction('&Follow Max Value', self)
        analyseMaxLineAction.triggered.connect(self.mainWidget._followMaxLineGUI)
        AnalyseMenu.addAction( analyseMaxLineAction)

        histoAction = qt.QAction('&Compute Histogram', self)
        histoAction.triggered.connect( self.mainWidget._histoGUI)
        AnalyseMenu.addAction(histoAction)

        fftcorrelationAction = qt.QAction('&Cross FFT Correlation', self)
        fftcorrelationAction.triggered.connect( self.mainWidget._FFTCORGUI)
        AnalyseMenu.addAction(fftcorrelationAction)

        mutualInformationAction = qt.QAction('&Mutual Information In Between Stack', self)
        mutualInformationAction.triggered.connect(self.mainWidget._MIGUI)
        AnalyseMenu.addAction(mutualInformationAction)

        ratioVolMaskAction = qt.QAction('&Local Volume Ratio Mask', self)
        ratioVolMaskAction.triggered.connect(self.mainWidget._RatioVolMaskGUI)
        AnalyseMenu.addAction(ratioVolMaskAction )

        ratioMaskAction = qt.QAction('&Local Ratio Mask 2D', self)
        ratioMaskAction.triggered.connect( self.mainWidget._RatioMaskGUI)
        AnalyseMenu.addAction(ratioMaskAction )

        SelectionMenu = ProcessMenu.addMenu('&Selection')
        smoothAction = qt.QAction('&Smooth Contours', self)
        smoothAction.triggered.connect( self.mainWidget._smoothContours)
        SelectionMenu.addAction(smoothAction)

        interpolateAction = qt.QAction('&Interpolate Contours', self)
        interpolateAction.triggered.connect(self.mainWidget._interpolateContourGUI)
        SelectionMenu.addAction(interpolateAction)

        interpolateMAction = qt.QAction('&Interpolate Mask', self)
        interpolateMAction.triggered.connect(self.mainWidget._interpolateMaskGUI)
        SelectionMenu.addAction(interpolateMAction)

        fftcorrelationAction = qt.QAction('&Cross FFT Correlation', self)
        fftcorrelationAction.triggered.connect( self.mainWidget._FFTCORGUI)
        AnalyseMenu.addAction(fftcorrelationAction)

        SegmentationMenu = ProcessMenu.addMenu('&Segmentation')

        segment_contour_Action = qt.QAction('&Segmentation From Contour', self)
        segment_contour_Action.triggered.connect( self.mainWidget._segFromContour)
        SegmentationMenu.addAction(segment_contour_Action)

        segment_th_Action = qt.QAction('&Threshold ', self)
        segment_th_Action.setStatusTip('3D Region Growing')
        segment_th_Action.triggered.connect( self.mainWidget._th_segmentation_GUI)
        SegmentationMenu.addAction(segment_th_Action)

        segment_rg_Action = qt.QAction('&3D Region Growing Threathold', self)
        segment_rg_Action.setStatusTip('3D Region Growing')
        segment_rg_Action.triggered.connect( self.mainWidget._rg_segmentation_GUI)
        SegmentationMenu.addAction(segment_rg_Action)

        segment_rgc_Action = qt.QAction('&3D Region Growing Confidence', self)
        segment_rgc_Action.setStatusTip('3D Region Growing')
        segment_rgc_Action.triggered.connect(self.mainWidget._rgc_segmentation_GUI)
        SegmentationMenu.addAction(segment_rgc_Action)

        segment_fm_Action = qt.QAction('&Fast Marching', self)
        segment_fm_Action.setStatusTip('Fast Marching Segmentation')
        segment_fm_Action.triggered.connect( self.mainWidget._fm_segmentation_GUI)
        SegmentationMenu.addAction(segment_fm_Action)

        segment_sd_Action = qt.QAction('&Shape Detection Level Set ITK', self)
        segment_sd_Action.setStatusTip('Shape Detection Level Set ITK')
        segment_sd_Action.triggered.connect(self.mainWidget._sd_segmentation_GUI)
        SegmentationMenu.addAction(segment_sd_Action)

        segment_geo_Action = qt.QAction('&Geodesic  Level Set ITK', self)
        segment_geo_Action.setStatusTip('Geodesic Level Set ITK')
        segment_geo_Action.triggered.connect(self.mainWidget._geo_segmentation_GUI)
        SegmentationMenu.addAction(segment_geo_Action)

        segment_snake_cv_Action = qt.QAction('&Morpho Snake Chan Vese', self)
        segment_snake_cv_Action .setStatusTip('Morpho Snake Chan Vese')
        segment_snake_cv_Action .triggered.connect(self.mainWidget._snake_chan_vese_GUI)
        SegmentationMenu.addAction(segment_snake_cv_Action )


        segment_wh_Action = qt.QAction('&Watershed', self)
        segment_wh_Action.setStatusTip('Watershed')
        segment_wh_Action.triggered.connect( self.mainWidget._wh_segmentation_GUI)
        SegmentationMenu.addAction(segment_wh_Action)

        EdgeDetectionMenu = ProcessMenu.addMenu('&Feature Detection')

        edge_GradientRGau_Action = qt.QAction('&Gradient Recursive Gaussian ', self)
        edge_GradientRGau_Action.setStatusTip('Gradient Recursive Gaussian ')
        edge_GradientRGau_Action.triggered.connect( self.mainWidget.edge_GradGauss_GUI)
        EdgeDetectionMenu.addAction(edge_GradientRGau_Action)

        edge_Canny_Action = qt.QAction('&Canny Edge', self)
        edge_Canny_Action.setStatusTip('Canny Edge')
        edge_Canny_Action.triggered.connect(self.mainWidget.edge_Canny_GUI)
        EdgeDetectionMenu.addAction(edge_Canny_Action)

        Zero_Crossing_Action = qt.QAction('&Zero Crossing Edge', self)
        Zero_Crossing_Action.setStatusTip('Zero Crossing Edge')
        Zero_Crossing_Action.triggered.connect( self.mainWidget.zero_Crossing_GUI)
        EdgeDetectionMenu.addAction(Zero_Crossing_Action)

        hessian_Action = qt.QAction('&Hessian Detection', self)
        hessian_Action.setStatusTip('Hessian Detection')
        hessian_Action.triggered.connect( self.mainWidget.hessian_GUI)
        EdgeDetectionMenu.addAction(hessian_Action)

        frangi_Action = qt.QAction('&Frangi Detection', self)
        frangi_Action.setStatusTip('Frangi Detection')
        frangi_Action.triggered.connect(self.mainWidget.frangi_GUI)
        EdgeDetectionMenu.addAction(frangi_Action)
 
        FilterMenu = ProcessMenu.addMenu('&Filter')

        filter_ani_diff_Action = qt.QAction('&Anisotropic Diffusion', self)
        filter_ani_diff_Action.setStatusTip('Anisotropic Diffusion')
        filter_ani_diff_Action.triggered.connect( self.mainWidget.filter_anidiff_GUI)
        FilterMenu.addAction(filter_ani_diff_Action)

        filter_recursiveGauss_Action = qt.QAction('&Recursive Gaussian', self)
        filter_recursiveGauss_Action.setStatusTip('Recursive Gaussian')
        filter_recursiveGauss_Action.triggered.connect( self.mainWidget.filter_recursiveGauss_GUI)
        FilterMenu.addAction(filter_recursiveGauss_Action)

        filter_median_Action = qt.QAction('&Median ', self)
        filter_median_Action.setStatusTip('Median Filter')
        filter_median_Action.triggered.connect( self.mainWidget.filter_median_GUI)
        FilterMenu.addAction(filter_median_Action)

        filter_sigmo_Action = qt.QAction('&Sigmoid Contrast ', self)
        filter_sigmo_Action.setStatusTip('Sigmoid Mapping  ')
        filter_sigmo_Action.triggered.connect(self.mainWidget.filter_Sigmo_GUI)
        FilterMenu.addAction(filter_sigmo_Action)


        filter_wv_Action = qt.QAction('&Wavelet Filtering', self)
        filter_wv_Action.setStatusTip('Wavelet Filtering')
        filter_wv_Action.triggered.connect( self.mainWidget.filter_WL_GUI)
        FilterMenu.addAction(filter_wv_Action)


        filter_mip_Action = qt.QAction('&Maximun/Minimum Image Projection', self)
        filter_mip_Action.setStatusTip('Maximum/Minimum Image Projection')
        filter_mip_Action.triggered.connect( self.mainWidget.filter_MIP_GUI)
        FilterMenu.addAction(filter_mip_Action)



        mathAction = qt.QAction( '&Math', self)
        mathAction.setStatusTip('Volume Simple Math Calculation')
        mathAction.triggered.connect(self.mainWidget._math_GUI)
        ProcessMenu.addAction(mathAction)

        ImageRegisteringMenu = ProcessMenu.addMenu('&Image Registration')
        registrationAction = qt.QAction(qt.QIcon('morpho.png'), '&Elastic Registration', self)
        registrationAction.triggered.connect( self.mainWidget._registration_GUI)
        ImageRegisteringMenu.addAction(registrationAction)

        morphoAction = qt.QAction(qt.QIcon('morpho.png'), '&Morphologic Operation', self)
        morphoAction.setStatusTip('Volume Morph Operation')
        morphoAction.triggered.connect(self.mainWidget._morpho_GUI)
        ProcessMenu.addAction(morphoAction)

        viewMenu = menubar.addMenu('&View')
        alphaAction = qt.QAction('&Alpha Map', self)
        alphaAction.setStatusTip('Add Colormap Overlay')
        alphaAction.triggered.connect(self.mainWidget.add_alpha_map_GUI)
        viewMenu.addAction(alphaAction)

        rmAlphaAction = qt.QAction('&Remove Overlay', self)
        rmAlphaAction.setStatusTip('Remove Overlay')
        rmAlphaAction.triggered.connect( self.mainWidget.remove_alpha_map)
        viewMenu.addAction(rmAlphaAction)


        StudyMenu = menubar.addMenu('&Study')
        venti_Action = qt.QAction('&Ventilation Map', self)
        venti_Action.setStatusTip('Compute Ventilation Map')
        venti_Action.triggered.connect(self.mainWidget.study_venti_GUI)
        StudyMenu.addAction(venti_Action)

        sVair_Action = qt.QAction('&Test sVair', self)
        sVair_Action.setStatusTip('Test sVair')
        sVair_Action.triggered.connect( self.mainWidget.sVair_GUI)
        StudyMenu.addAction(sVair_Action)

        density_Action = qt.QAction('&Density Change', self)
        density_Action.setStatusTip('Compute Density Change')
        density_Action.triggered.connect(self.mainWidget.study_density_GUI)
        StudyMenu.addAction(density_Action)

if __name__ == "__main__":


    app = qt.QApplication(["-display"])
    foo = MyMainWindow()
    foo.resize(400, 400)
    foo.show()
    sys.exit(app.exec_())
