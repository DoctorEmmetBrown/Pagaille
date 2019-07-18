# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 13:30:03 2017

@author: BROCHE

Modified by FARDIN. Added Saving and Chess Image Control.
"""
import SimpleITK as sitk
import numpy as np
import glob
import EdfFile
import os


def change_grid_size(R):
    print('New resolution')


class Registering():

    def __init__(self):
        print('Init Registration')
        self.reg_method = sitk.ImageRegistrationMethod()
        # self.reg_method.SetNumberOfThreads(32)

    def load_image(self, array_fixed, array_moving):
        print('Init Image')
        self.i_fixed = sitk.GetImageFromArray(array_fixed, isVector=False)
        self.i_moving = sitk.GetImageFromArray(array_moving, isVector=False)

    def init_bspline(self, gridSize, order):
        print('Init Bspline')
        x_grid_size = gridSize[0]
        y_grid_size = gridSize[1]
        z_grid_size = gridSize[2]

        grid_physical_spacing = [x_grid_size, y_grid_size, z_grid_size]
        image_physical_size = [size * spacing for size, spacing in
                               zip(self.i_fixed.GetSize(), self.i_fixed.GetSpacing())]

        mesh_size = [int(image_size / grid_spacing + 0.5) for image_size, grid_spacing in
                     zip(image_physical_size, grid_physical_spacing)]
        self.initial_transform = sitk.BSplineTransformInitializer(image1=self.i_fixed,
                                                                  transformDomainMeshSize=mesh_size, order=order)
        self.reg_method.SetInitialTransform(self.initial_transform)

    def init_affine(self):
        self.affine = sitk.AffineTransform(self.i_fixed.GetDimension())
        self.affine.Scale(1 / 0.9, True)
        self.reg_method.SetInitialTransform(self.affine)

    def init_Translation(self):
        self.translate = sitk.TranslationTransform(self.i_fixed.GetDimension())
        self.reg_method.SetInitialTransform(self.translate)

    def init_Euler(self):
        self.Euler = sitk.Euler3DTransform()
        self.reg_method.SetInitialTransform(self.Euler)

    def SetInitialTransform(self, T):

        T1 = self.reg_method.GetInitialTransform()
        T2 = sitk.Transform(T1)
        T2.AddTransform(T)
        self.reg_method.SetInitialTransform(T2)

    def SetScalableMesh(self, InitialMesh, factors):
        self.reg_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, lambda: change_grid_size(self.reg_method))

    def init_metric(self, Metric, parameters=[], BoolGradient=False, SamplingPerc=0.5):
        print('Init Metric')
        if Metric == "Means Squares":
            self.reg_method.SetMetricAsMeanSquares()

        elif Metric == "Correlation":
            self.reg_method.SetMetricAsCorrelation()

        elif Metric == "Joint Histogram Mutual Information":
            self.reg_method.SetMetricAsJointHistogramMutualInformation(parameters[0], parameters[1])

        elif Metric == "Mattes Mutual Information":
            self.reg_method.SetMetricAsMattesMutualInformation(parameters[0])

        elif Metric == "Neighborhood Correlation (ANTs)":
            self.reg_method.SetMetricAsANTSNeighborhoodCorrelation(parameters[0])

        self.reg_method.SetMetricSamplingStrategy(self.reg_method.RANDOM)
        self.reg_method.SetMetricSamplingPercentage(SamplingPerc)

        self.reg_method.SetMetricUseFixedImageGradientFilter(BoolGradient)
        self.reg_method.SetMetricUseMovingImageGradientFilter(BoolGradient)

    def initOptimizer(self, Optimizer, Parameters):
        print('Init Optimizer')
        """ Optimizer """
        if Optimizer == "Regular Step Gradient Descent":
            par1 = Parameters[0]
            par2 = Parameters[1]
            par3 = Parameters[2]
            par4 = Parameters[3]
            par5 = Parameters[4]

            if Parameters[5] == 0:
                par6 = self.reg_method.Never
            elif Parameters[5] == 1:
                par6 = self.reg_method.Once
            elif Parameters[5] == 2:
                par6 = self.reg_method.EachIteration

            par7 = Parameters[6]
            self.reg_method.SetOptimizerAsRegularStepGradientDescent(par1, par2, par3, par4, par5, par6, par7)

        elif Optimizer == "Gradient Descent":
            par1 = Parameters[0]
            par2 = Parameters[1]
            par3 = Parameters[2]
            par4 = Parameters[3]

            if Parameters[4] == 0:
                par5 = self.reg_method.Never
            elif Parameters[4] == 1:
                par5 = self.reg_method.Once
            elif Parameters[4] == 2:
                par5 = self.reg_method.EachIteration
            par6 = Parameters[5]
            self.reg_method.SetOptimizerAsGradientDescent(par1, par2, par3, par4, par5, par6)

        elif Optimizer == "Gradient Descent Line Search":

            par1 = float(Parameters[0])
            par2 = int(Parameters[1])
            par3 = float(Parameters[2])
            par4 = int(Parameters[3])
            par5 = float(Parameters[4])
            par6 = float(Parameters[5])
            par7 = float(Parameters[6])
            par8 = int(Parameters[7])
            if int(Parameters[8]) == 0:
                par9 = self.reg_method.Never
            elif int(Parameters[8]) == 1:
                par9 = self.reg_method.Once
            elif int(Parameters[8]) == 2:
                par9 = self.reg_method.EachIteration
            par10 = int(Parameters[9])
            self.reg_method.SetOptimizerAsGradientDescentLineSearch(par1, par2, par3, par4, par5, par6, par7, par8,
                                                                    par9, par10)

        elif Optimizer == "Conjugate Gradient Line Search":
            par1 = float(Parameters[0])
            par2 = int(Parameters[1])
            par3 = float(Parameters[2])
            par4 = int(Parameters[3])
            par5 = float(Parameters[4])
            par6 = float(Parameters[5])
            par7 = float(Parameters[6])
            par8 = int(Parameters[7])
            if int(Parameters[8]) == 0:
                par9 = self.reg_method.Never
            elif int(Parameters[8]) == 1:
                par9 = self.reg_method.Once
            elif int(Parameters[8]) == 2:
                par9 = self.reg_method.EachIteration
            par10 = int(Parameters[9])
            self.reg_method.SetOptimizerAsConjugateGradientLineSearch(par1, par2, par3, par4, par5, par6, par7, par8,
                                                                      par9, par10)

        elif Optimizer == "LBFGSB":
            par1 = Parameters[0]
            par2 = Parameters[1]
            par3 = Parameters[2]
            par4 = Parameters[3]
            par5 = Parameters[4]
            self.reg_method.SetOptimizerAsLBFGSB(par1, par2, par3, par4, par5)

        elif Optimizer == "Powell":
            par1 = Parameters[0]
            par2 = Parameters[1]
            par3 = Parameters[2]
            par4 = Parameters[3]
            par5 = Parameters[4]
            self.reg_method.SetOptimizerAsPowell(par1, par2, par3, par4, par5)

        elif Optimizer == "Amoeba":
            par1 = Parameters[0]
            par2 = Parameters[1]
            par3 = Parameters[2]
            par4 = Parameters[3]
            self.reg_method.SetOptimizerAsAmoeba(par1, par2, par3, par4)

    def init_Interpolator(self, Interpolator):
        print('Init Interpolator')
        if Interpolator == "Nearest neighbor":
            self.interpolator = sitk.sitkNearestNeighbor
        elif Interpolator == "Linear Interpolation":
            self.interpolator = sitk.sitkLinear
        elif Interpolator == "BSpline":
            self.interpolator = sitk.sitkBSpline
        elif Interpolator == "Gaussian":
            self.interpolator = sitk.sitkGaussian
        elif Interpolator == "Label Gaussian":
            self.interpolator = sitk.sitkLabelGaussian
        elif Interpolator == "Hamming Windowed Sinc":
            self.interpolator = sitk.sitkHammingWindowedSinc
        elif Interpolator == "Cosine Windowed Sinc":
            self.interpolator = sitk.sitkCosineWindowedSinc
        elif Interpolator == "Welch Windowed Sinc":
            self.interpolator = sitk.sitkWelchWindowedSinc
        elif Interpolator == "Lanczos Windowed Sinc":
            self.interpolator = sitk.sitkLanczosWindowedSinc
        elif Interpolator == "Blackman Windowed Sinc":
            self.interpolator = sitk.sitkBlackmanWindowedSinc

        self.reg_method.SetInterpolator(self.interpolator)

    def initScaling(self, scalingVector, smoothingVector):
        print('Init Scaling')
        self.reg_method.SetShrinkFactorsPerLevel(scalingVector)
        self.reg_method.SetSmoothingSigmasPerLevel(smoothingVector)
        self.reg_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    def initMask(self, InputMaskF, InputMaskM):
        print("Init Masks")
        self.MF = sitk.GetImageFromArray(InputMaskF, isVector=False)
        self.MM = sitk.GetImageFromArray(InputMaskM, isVector=False)
        self.reg_method.SetMetricFixedMask(self.MF)
        self.reg_method.SetMetricMovingMask(self.MM)

    def SetZeroTransform(self, T0):

        # self.Transform = sitk.DisplacementFieldTransform(T0.GetDimension()) #self.i_transform
        # self.Transform.SetDisplacementField(T0)

        self.reg_method.SetMovingInitialTransform(T0)

    def Execute(self):
        print('Runing Registration ...')
        print('Initial value :')
        print(self.reg_method.MetricEvaluate(self.i_fixed, self.i_moving))
        self.Final_Transform = self.reg_method.Execute(self.i_fixed, self.i_moving)
        print(self.reg_method.GetOptimizerStopConditionDescription())

    def GetTransform(self):
        print(self.Final_Transform)
        return self.Final_Transform

    def computeJacobian(self):
        print('Compute Jacobian')
        filterJacob = sitk.DisplacementFieldJacobianDeterminantFilter()
        return sitk.ArrayFromImage(filterJacob.Execute(self.displacementField))

    def applyFinalTransform(self, ImageIn):
        print('Apply Transform to external Image')
        ImageITK = sitk.GetImageFromArray(ImageIn)
        ImageOut = sitk.Resample(ImageITK, self.ImageFixe, self.interpolator, 0.0, ImageITK.GetPixelIDValue())
        return np.copy(sitk.GetArrayFromImage(ImageOut))

    def Apply_Step(self, image):
        print('Apply Transform to Image Moving')
        i_moving = sitk.GetImageFromArray(image, isVector=False)
        i_moved = sitk.Resample(i_moving, self.i_fixed, self.Final_Transform, self.interpolator, 0.0,
                                self.i_moving.GetPixelIDValue())
        return np.copy(sitk.GetArrayFromImage(i_moved))

    def ApplyMovingTransform(self):
        print('Apply Transform to Image Moving')
        self.i_fixed = sitk.GetImageFromArray(image_fixed)
        self.i_moved = sitk.Resample(self.i_moving, self.i_fixed, self.Final_Transform, self.interpolator, 0.0,
                                     self.i_moving.GetPixelIDValue())

    def computeVectorField(self):
        print('Coverting Final Bspline to Vector Field')
        FilterTransform = sitk.TransformToDisplacementFieldFilter()
        FilterTransform.SetReferenceImage(self.i_fixed)
        self.displacementField = FilterTransform.Execute(self.Final_Transform)

        return sitk.ArrayFromImage(self.displacementField)

    def SaveTransform(self, output_Transform):

        # Check Folder Exists
        if not os.path.exists(output_Transform):
            os.makedirs(output_Transform)

        # Write Transform
        File_Transform = output_Transform + "Transform.tfm"
        sitk.WriteTransform(self.Final_Transform, File_Transform)
        print(self.reg_method.GetMetricValue())

    def SaveCompositeTransform(self, output_Transform):

        # Check Folder Exists
        if not os.path.exists(output_Transform):
            os.makedirs(output_Transform)

        # Write Transform
        File_Transform = output_Transform + "Transform.tfm"

        sitk.WriteTransform(self.Transform.AddTransform(self.Final_Transform), File_Transform)
        print(self.reg_method.GetMetricValue())


if __name__ == "__main__":

    path_fixed = '/data/id17/broncho/IH2887/Reconstruction/rabbit2_sc6_/10resized/'
    path_moving = '/data/id17/broncho/IH2887/Reconstruction/rabbit2_sc6_/200resized/'
    path_maskf = '/data/id17/broncho/IH2887/Reconstruction/rabbit2_sc6_/10mask/'
    path_maskm = '/data/id17/broncho/IH2887/Reconstruction/rabbit2_sc6_/200mask/'
    output_Transform = '/data/id17/broncho/IH2887/Reconstruction/rabbit2_sc6_/10_250_transform/'

    print("Loading Images ...")
    image_fixed_names = glob.glob(path_fixed + "/*.edf")
    image_fixed_names = sorted(image_fixed_names)
    fileEdf = EdfFile.EdfFile(image_fixed_names[0], access='rb')
    image = fileEdf.GetData(0)
    image_fixed = np.zeros((image.shape[0], image.shape[1], len(image_fixed_names)))

    z = 0
    for name_image in image_fixed_names:
        print(name_image)
        fileEdf = EdfFile.EdfFile(name_image, access='rb')
        image_fixed[:, :, z] = fileEdf.GetData(0)
        z += 1

    image_moving_names = glob.glob(path_moving + "/*.edf")
    image_moving_names = sorted(image_moving_names)
    fileEdf = EdfFile.EdfFile(image_moving_names[0], access='rb')
    image = fileEdf.GetData(0)
    image_moving = np.zeros((image.shape[0], image.shape[1], len(image_moving_names)))

    z = 0
    for name_image in image_moving_names:
        print(name_image)
        fileEdf = EdfFile.EdfFile(name_image, access='rb')
        image_moving[:, :, z] = fileEdf.GetData(0)
        z += 1

    mask_fixed_names = glob.glob(path_maskf + "/*.edf")
    mask_fixed_names = sorted(mask_fixed_names)
    fileEdf = EdfFile.EdfFile(mask_fixed_names[0], access='rb')
    image = fileEdf.GetData(0)
    InputMaskF = np.zeros((image.shape[0], image.shape[1], len(mask_fixed_names)))

    z = 0
    for name_image in mask_fixed_names:
        print(name_image)
        fileEdf = EdfFile.EdfFile(name_image, access='rb')
        InputMaskF[:, :, z] = fileEdf.GetData(0)
        z += 1

    mask_moving_names = glob.glob(path_maskm + "/*.edf")
    mask_moving_names = sorted(mask_moving_names)
    fileEdf = EdfFile.EdfFile(mask_moving_names[0], access='rb')
    image = fileEdf.GetData(0)
    InputMaskM = np.zeros((image.shape[0], image.shape[1], len(mask_moving_names)))

    z = 0
    for name_image in mask_moving_names:
        print(name_image)
        fileEdf = EdfFile.EdfFile(name_image, access='rb')
        InputMaskM[:, :, z] = fileEdf.GetData(0)
        z += 1

    reg_method = Registering()
    reg_method.load_image(image_fixed, image_moving)
    reg_method.initMask(InputMaskF, InputMaskM)
    reg_method.init_bspline([40, 40, 40], 3)
    reg_method.init_metric("Means Squares", [], True, 0.5)
    reg_method.init_Interpolator("BSpline")
    reg_method.initOptimizer("Conjugate Gradient Line Search", [10.0, 1000, 1e-06, 10, 0.01, 10.0, 0.01, 20.0, 2, 10])
    reg_method.initScaling([16, 16, 8, 4], [0, 0, 0, 0])
    reg_method.Execute()
    reg_method.SaveTransform(output_Transform)
    # vector_field = reg_method.computeVectorField()
    # jacobian = reg_method.computeJacobian()

