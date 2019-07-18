'''

Based on ImageProcessing.py in PINI3 by L. Broche

'''

import InputOutput_fabio as InputOutput
import SimpleITK as sitk
import math
import numpy as np
import sys


# import ImageProcessing


class MaskCreator():

    def __init__(self):
        print('Init Mask Creator')
        self.IO = InputOutput.InputOutput()

    def LoadImages(self, path):

        image_in = self.IO.ReadEdfVolume(path)
        ITK_Vol = sitk.GetImageFromArray(image_in, isVector=False)
        return ITK_Vol

    def SaveImages(self, ITKVolume, path, root):

        print('Saving')

        ImageToSave = sitk.GetArrayFromImage(ITKVolume)
        self.IO.SaveEdf(ImageToSave, path, root)

    def SaveEdf(self, ImageToSave, path, root):

        print('Saving')
        self.IO.SaveEdf(ImageToSave, path, root)

    def resizeImage(self, img, interpolator, spacing):

        print('Resize')

        if interpolator == "Nearest neighbor":
            interpolatorITK = sitk.sitkNearestNeighbor
        elif interpolator == "Linear":
            interpolatorITK = sitk.sitkLinear
        elif interpolator == "BSpline":
            interpolatorITK = sitk.sitkBSpline
        elif interpolator == "Gaussian":
            interpolatorITK = sitk.sitkGaussia
        elif interpolator == "Label Gaussian":
            interpolatorITK = sitk.sitkLabelGaussian
        elif interpolator == "Hamming Windowed Sinc":
            interpolatorITK = sitk.sitkHammingWindowedSinc
        elif interpolator == "Cosine Windowed Sinc":
            interpolatorITK = sitk.sitkCosineWindowedSinc
        elif interpolator == "Welch Windowed Sinc":
            interpolatorITK = sitk.sitkWelchWindowedSinc
        elif interpolator == "Lanczos Windowed Sinc":
            interpolatorITK = sitk.sitkLanczosWindowedSinc
        elif interpolator == "Blackman Windowed Sinc":
            interpolatorITK = sitk.sitkBlackmanWindowedSinc

        if len(spacing) != img.GetDimension(): raise Exception("len(spacing) != " + str(img.GetDimension()))

        inSpacing = img.GetSpacing()
        inSize = img.GetSize()
        size = [int(math.ceil(inSize[i] * (inSpacing[i] / spacing[i]))) for i in range(img.GetDimension())]

        identityTransform = sitk.Transform()
        img = sitk.Resample(img, size, identityTransform, interpolatorITK, [0] * 3, spacing)

        return img

    def anisotropic_diffusion(self, ITK_Vol, time_step, conductance, nbIter):

        print('Filtering')
        curvDiff = sitk.CurvatureAnisotropicDiffusionImageFilter()
        curvDiff.SetTimeStep(time_step)
        curvDiff.SetConductanceParameter(conductance)
        curvDiff.SetNumberOfIterations(nbIter)
        ITK_Vol = curvDiff.Execute(ITK_Vol)
        return ITK_Vol

    def Seeds(self, filtered_img, segmentation_values, roi):
        filtered = sitk.GetArrayFromImage(filtered_img)
        print('Placing Seeds')
        seeds_x = np.random.randint(roi[0], roi[0] + roi[2], size=20)
        seeds_y = np.random.randint(roi[1], roi[1] + roi[3], size=20)
        seed_list = []
        k = 0
        for i in range(20):
            x = np.intp(seeds_x[i])
            y = np.intp(seeds_y[i])
            value = filtered[y, x, 0]
            if (value > segmentation_values[0] and value < segmentation_values[1]):
                seed_list.append([int(seeds_x[i]), int(seeds_y[i]), 0])
                k += 1
        print('Placed ' + str(k) + ' seeds')
        return seed_list

    def SegConnectedThreshold(self, ITK_Vol, seedListToSegment, segmentation_values):
        print('Segmenting')
        val_min = segmentation_values[0]
        val_max = segmentation_values[1]
        segmentationFilter = sitk.ConnectedThresholdImageFilter()

        for seed in seedListToSegment:
            seedItk = ([seed[2], seed[0], seed[1]])
            segmentationFilter.AddSeed(seedItk)

        segmentationFilter.SetLower(val_min)
        segmentationFilter.SetUpper(val_max)
        segmentationFilter.SetReplaceValue(1)
        ITK_Vol = segmentationFilter.Execute(ITK_Vol)
        image = sitk.GetArrayFromImage(ITK_Vol)

        number_of_pixel = np.sum(image, axis=None)

        return image.astype(np.uint8)

    def Close(self, vol, kernel):
        print('Closing')
        vol_itk = sitk.GetImageFromArray(vol, isVector=False)
        filterMorpho = sitk.BinaryMorphologicalClosingImageFilter()
        filterMorpho.SetKernelRadius(kernel)
        vol_out_itk = filterMorpho.Execute(vol_itk)
        return sitk.GetArrayFromImage(vol_out_itk)

    def Dilatation(self, vol, kernel):
        print('Dilatation')
        vol_itk = sitk.GetImageFromArray(vol, isVector=False)
        filterMorpho = sitk.BinaryDilateImageFilter()
        filterMorpho.SetKernelRadius(kernel)
        vol_out_itk = filterMorpho.Execute(vol_itk)
        return sitk.GetArrayFromImage(vol_out_itk)

    def Erosion(self, vol, kernel):
        print('Erosion')
        vol_itk = sitk.GetImageFromArray(vol, isVector=False)
        filterMorpho = sitk.BinaryErodeImageFilter()
        filterMorpho.SetKernelRadius(kernel)
        vol_out_itk = filterMorpho.Execute(vol_itk)
        return sitk.GetArrayFromImage(vol_out_itk)

    def Fill(self, vol):
        print('Filling')
        shapeVol = np.shape(vol)
        for z in range(0, shapeVol[2]):
            image = vol[:, :, z]
            imageITK = sitk.GetImageFromArray(image, isVector=False)
            filterMorpho = sitk.BinaryFillholeImageFilter()
            imageITK = filterMorpho.Execute(imageITK)
            image = sitk.GetArrayFromImage(imageITK)
            vol[:, :, z] = image
        return vol.astype(np.uint8)


if __name__ == "__main__":

    rabbit = '7'
    scan = '21'
    # delay_arr=[10,50,75,100,125,150,175,200,225,250,275,300,325,350]
    delay_arr = [10, 150]
    phase = '3'
    experiment = 'md1072'
    for delay in delay_arr:
        path = '/data/id17/broncho/' + experiment + '/Reconstruction/rabbit' + rabbit + '_sc' + scan + '/' + str(
            delay).zfill(3) + '_' + phase + '/'
        # path =  '/data/id17/broncho/IH2779/r2_sc_20_/Slices/'+str(delay).zfill(3)+'_07/'
        # path =  '/data/id17/broncho/md1072/Reconstruction/rabbit'+rabbit+'_sc'+scan+'/'+str(delay).zfill(3)+'_'+phase
        # path =  '/data/id17/broncho/md1072/Reconstruction/rabbit'+rabbit+'_sc'+scan+'/'+str(delay).zfill(3)
        # resized_path = '/data/id17/broncho/IH2887/Reconstruction/rabbit'+rabbit+'_sc'+scan+'_/'+str(delay)+'resizedm/'
        # resized_path = '/data/id17/broncho/IH2779/r2_sc_20_/Slices/'+str(delay)+'resized/'
        # resized_path = '/data/id17/broncho/md1072/Reconstruction/rabbit'+rabbit+'_sc'+scan+'/'+str(delay)+'resized/'
        # filtered_path = '/data/id17/broncho/IH2887/Reconstruction/rabbit'+rabbit+'_sc'+scan+'_/'+str(delay)+'filtered/'
        # filtered_path = '/data/id17/broncho/IH2779/r2_sc_20_/Slices/'+str(delay)+'filtered/'
        filtered_path = '/data/id17/broncho/' + experiment + '/Reconstruction/rabbit' + rabbit + '_sc' + scan + '/' + str(
            delay) + '_' + phase + 'filtered/'
        # filtered_path = '/data/id17/broncho/md1072/Reconstruction/rabbit'+rabbit+'_sc'+scan+'/'+str(delay)+'filtered/'
        # resized_root = 'rabbit'+rabbit+'_sc'+scan+'_'+str(delay).zfill(3)+'_'+phase+'_'
        resized_root = 'rabbit' + rabbit + '_sc' + scan + '_' + str(delay).zfill(3) + '_' + phase + '_'
        # mask_path = '/data/id17/broncho/IH2887/Reconstruction/rabbit'+rabbit+'_sc'+scan+'_/'+str(delay)+'mask/'
        # mask_path = '/data/id17/broncho/IH2779/r2_sc_20_/Slices/'+str(delay)+'mask/'
        mask_path = '/data/id17/broncho/' + experiment + '/Reconstruction/rabbit' + rabbit + '_sc' + scan + '/' + str(
            delay) + '_' + phase + 'mask/'
        mask_root = 'rabbit' + rabbit + '_sc' + scan + '_' + str(delay) + '_' + phase + '_m_'
        # print(path,resized_path,filtered_path,resized_root,mask_path,mask_root)

        mask_method = MaskCreator()

        in_image = mask_method.LoadImages(path)
        in_image = sitk.GetImageFromArray(sitk.GetArrayFromImage(in_image), isVector=False)
        res_image = mask_method.resizeImage(in_image, "BSpline", [2, 2, 2])
        # del in_image
        ##mask_method.SaveEdf(res_image,resized_path,resized_root)
        # mask_method.SaveImages(res_image,resized_path,resized_root)
        sys.stdout.flush()

        # res_image = mask_method.LoadImages(resized_path)
        #	res_image=sitk.GetImageFromArray(res_image,isVector=False)
        filt_image = mask_method.anisotropic_diffusion(res_image, 0.06, 9.0, 25)  # 0.06,9.0,50
        mask_method.SaveImages(filt_image, filtered_path, resized_root)

    # filt_image=mask_method.LoadImages(filtered_path)
    # print filt_image.GetSize()
# x=1833
# y=1141
# print filt_image[0,x,y]

# Now place some seeds in the image
# roi=[732,700,1108,1100]#[548,620,1254,1332]  #[732,700,1108,1100]
# segmentation_values=[-0.35,0.1]
# seeds=mask_method.Seeds(filt_image,segmentation_values,roi)

# Run segmentation

# mask_img=mask_method.SegConnectedThreshold(filt_image,seeds,segmentation_values)
# mask_img=mask_method.Close(mask_img,3)
##mask_img=mask_method.Dilatation(mask_img,5)
##mask_img=mask_method.Fill(mask_img)
# mask_method.SaveEdf(mask_img,mask_path,mask_root)

# filt_image=mask_method.LoadImages(filtered_path)
# mask_image=mask_method.Segment(filt_image,)

# med = sitk.MedianImageFilter()
# med.SetRadius(5)
# ITK_Vol = med.Execute(res_image)

# filtering=sitk.GradientMagnitudeImageFilter()
# output_image=filtering.Execute(ITK_Vol)
# mask_method.SaveImages(output_image,filtered_path,resized_root)

# [1200:1200+3200,1100:1100+3400]
