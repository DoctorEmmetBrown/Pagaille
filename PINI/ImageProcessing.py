# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:26:20 2016

@author: broche
"""

import numpy as np
import PyMcaQt as qt
import SimpleITK as sitk

import scipy
from scipy.optimize import curve_fit
from scipy import stats
from scipy.interpolate import splprep, splev

from scipy import ndimage
from scipy.ndimage.measurements import label


from skimage import measure

import math
import pywt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from time import gmtime, strftime
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

import entropy_estimators as ee
from numpy import linalg as LA

from itertools import combinations_with_replacement
from scipy import ndimage as ndi

def divide_nonzero(array1, array2):
    """
    Divides two arrays. Returns zero when dividing by zero.
    """
    denominator = np.copy(array2)
    denominator[denominator == 0] = 1e-10
    return np.divide(array1, denominator)


def create_image_like(data, image):
    return image.__class__(data, affine=image.affine, header=image.header)


def absolute_eigenvaluesh(nd_array):
    """
    Computes the eigenvalues sorted by absolute value from the symmetrical matrix.
    :param nd_array: array from which the eigenvalues will be calculated.
    :return: A list with the eigenvalues sorted in absolute ascending order (e.g. [eigenvalue1, eigenvalue2, ...])
    """
    eigenvalues = np.linalg.eigvalsh(nd_array)
    sorted_eigenvalues = sortbyabs(eigenvalues, axis=-1)
    return [np.squeeze(eigenvalue, axis=-1)
            for eigenvalue in np.split(sorted_eigenvalues, sorted_eigenvalues.shape[-1], axis=-1)]


def sortbyabs(a, axis=0):
    """Sort array along a given axis by the absolute value
    modified from: http://stackoverflow.com/a/11253931/4067734
    """
    index = list(np.ix_(*[np.arange(i) for i in a.shape]))
    index[axis] = np.abs(a).argsort(axis)
    return a[index]




def equalizeAlongLine(Image, posXY):
    Image = Image.astype(float)
    mean_array = []
    for pos in posXY:
        x = int(pos[0])
        y = int(pos[1])
        z = int(pos[2])

        roi = Image[z, x - 5:x + 5, y - 5:y + 5]
        mean_array.append(np.mean(np.mean(roi)))

    mean_value = np.mean(mean_array)

    for z in range(0, Image.shape[0]):
        Image[z, :, :] = np.add(Image[z, :, :], (mean_value - mean_array[z]))

    return Image


def computeVessels(image):
    labeled_array, num_features = label(image)

    props = measure.regionprops(labeled_array)

    areas = []

    for prop in props:
        circular = (np.pi * 4.0 * prop.area) / (prop.perimeter ** 2.0)
        if (circular > 0.8) and (circular < 1.2):
            areas.append(prop.area * (0.25 * 0.25))

    return areas


def BackSigmoid(Image, Alpha, Beta):
    Image[Image <= 0.0593] = 0.0593
    Image[Image >= 224.60] = 224.60


    Image = -Alpha * (np.log((224.61 / Image) - 1.0)) + Beta
    return Image


def compute_sVair(ImageFix, ImageMov, Trachea, Distance):
    Trachea[Trachea == 1] = -1
    Trachea += 1
    ImageFix *= Trachea
    ImageMov *= Trachea

    # Grad = gradGauss(ImageFix,0.2)

    # ImageFix = BackSigmoid(ImageFix,-100,-824)
    # ImageMov = BackSigmoid(ImageMov,-100,-824)
    ImageMov[ImageMov == ImageMov[0, 0, 0]] = -2000
    ImageFix[ImageFix == ImageFix[0, 0, 0]] = -2000

    # ImageMov = np.copy(grey_dilation(ImageMov, size=(40,40,40)))
    # ImageMov2[ImageMov > 2.0] = 0.0
    # ImageMov = ImageMov2 + ImageMov

    # ImageMov = gaussian_filter(ImageMov,5)

    # ImageMov = np.nan_to_num(ImageMov)
    # ImageMov[ImageMov >= 0] = -0.1
    # ImageMov[ImageMov <= -1000] = -999.9

    # ImageFix = np.nan_to_num(ImageFix)
    # ImageFix[ImageFix >= 0] = -0.1
    # ImageFix[ImageFix <= -1000] = -999.9

    # sVair = 1000.0*((ImageMov-ImageFix)/(ImageFix*(ImageMov+1000.0)))


    ImageMov = median(ImageMov, 1,False)
    ImageFix = median(ImageFix, 1, False)

    # ImageMov = grey_dilation(ImageMov, size=(5,5,5))

    deltaHU = (ImageMov - ImageFix)
    deltaHU[ImageFix < -1100] = 0.0
    deltaHU[ImageFix > 200] = 0.0
    #deltaHU[deltaHU < -100.0] = 0
    #deltaHU[deltaHU > 500.0] = 0

    sVair = -1000.0 * (deltaHU) / (ImageFix * (ImageMov + 1000.0))

    sVair[ImageFix < -1100] = 0.0
    sVair[ImageFix > 200] = 0.0
    #sVair[sVair < -1.0] = 0.0
    #sVair[sVair > 5.0] = 0.0

    Distance[ImageFix < -1100] = 0.0
    Distance[ImageFix > 200] = 0.0
    #Distance[deltaHU < -500.0] = 0.0
    #Distance[deltaHU > 500.0] = 0.0

    Cluster = np.zeros((ImageFix.shape[0], ImageFix.shape[1], ImageFix.shape[2]))

    Cluster[(ImageMov < -900) * (ImageFix < -900)] = 5
    Cluster[(ImageMov > -900) * (ImageMov < -500) * (ImageFix < -900)] = 1
    Cluster[(ImageMov > -500) * (ImageFix < -900)] = 6

    Cluster[(ImageMov < -900) * (ImageFix < -500) * (ImageFix > -900)] = 7
    Cluster[(ImageMov > -900) * (ImageMov < -500) * (ImageFix < -500) * (ImageFix > -900)] = 1
    Cluster[(ImageMov > -500) * (ImageFix < -500) * (ImageFix > -900)] = 9

    Cluster[(ImageMov < -900) * (ImageFix > -500)] = 3
    Cluster[(ImageMov > -900) * (ImageMov < -500) * (ImageFix > -500)] = 2
    Cluster[(ImageMov > -500) * (ImageFix > -500)] = 8


    #Cluster[(ImageMov > -832) * (ImageFix > -926)] = 2
    #Cluster[(ImageMov > -832) * (ImageFix < -926)] = 1
    #Cluster[(ImageMov < -832) * (ImageFix > -926)] = 3
    #Cluster[(ImageMov < -832) * (ImageFix < -926)] = 4
    #Cluster[(ImageMov > -530)] = 1
    #Cluster[(ImageFix > -530) ] = 5

    Cluster[ImageFix < -1100] = 0.0
    Cluster[ImageFix > 200] = 0.0

    Cluster[0, 0, 0] = 1
    Cluster[0, 0, 1] = 2
    Cluster[0, 1, 0] = 3
    Cluster[0, 1, 1] = 4
    Cluster[1, 0, 0] = 5
    Cluster[1, 0, 1] = 6
    Cluster[1, 1, 0] = 7
    Cluster[1, 1, 1] = 8
    Cluster[0, 0, 2] = 9

    # Distance -=  ((Grad/np.max(Grad))*np.max(Distance)*0.3)
    # sVair -= ((Grad/np.max(Grad))*np.max(sVair)*0.3)
    # deltaHU -= ((Grad/np.max(Grad))*np.max(deltaHU)*0.3)
    # ImageMov -= ((Grad/np.max(Grad))*np.max(ImageMov)*0.3)
    # ImageFix -= ((Grad/np.max(Grad))*np.max(ImageFix)*0.3)

    return [deltaHU, sVair, Distance, Cluster, ImageMov, ImageFix]


def mip(inputDataToSeg, mip, direction, thickness):
    if direction == 'x':
        ratio = int(inputDataToSeg.shape[0] / thickness)
        image_out = np.zeros((ratio, inputDataToSeg.shape[1], inputDataToSeg.shape[2]))
    elif direction == 'y':
        ratio = int(inputDataToSeg.shape[1] / thickness)
        image_out = np.zeros((inputDataToSeg.shape[0], ratio, inputDataToSeg.shape[2]))
    elif direction == 'z':
        ratio = int(inputDataToSeg.shape[2] / thickness)
        image_out = np.zeros((inputDataToSeg.shape[0], inputDataToSeg.shape[1], ratio))

    for i in range(0, ratio-1):
        if direction == 'x':
            if mip == 'Min':
                image_out[i, :, :] = np.min(inputDataToSeg[i * thickness:(i + 1) * thickness, :, :], axis=0)
            elif mip == 'Max':
                image_out[i, :, :] = np.max(inputDataToSeg[i * thickness:(i + 1) * thickness, :, :], axis=0)
        elif direction == 'y':
            if mip == 'Min':
                image_out[:, i, :] = np.min(inputDataToSeg[:, i * thickness:(i + 1) * thickness, :], axis=1)
            elif mip == 'Max':
                image_out[:, i, :] = np.max(inputDataToSeg[:, i * thickness:(i + 1) * thickness, :], axis=1)
        elif direction == 'z':
            if mip == 'Min':
                image_out[:, :, i] = np.min(inputDataToSeg[:, :, i * thickness:(i + 1) * thickness], axis=2)
            elif mip == 'Max':
                image_out[:, :, i] = np.max(inputDataToSeg[:, :, i * thickness:(i + 1) * thickness], axis=2)

    return image_out


def gauss(x, mu, sigma, A):
    return A * math.exp(-(x - mu) ** 2 / 2 / sigma ** 2)


def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)


def equalizeHisto(vol):
    data = vol[vol.shape[0] / 2, :, :]
    sizeData = data.shape[0] * data.shape[1]
    data = np.random.choice(np.ravel(data), sizeData / 10)
    data = data / float(data.max())

    range = [data.min(), data.max()]
    bins = np.arange(range[0], range[1] + 2, (range[1] - range[0]) / 1000.0) - 0.5
    heights, edges = np.histogram(data, bins=bins, range=range)
    edges = edges[:-1] + (edges[1] - edges[0])

    expected = (np.mean(data) * 1.2, np.var(data), sizeData / 2400, np.mean(data) * 0.8, np.var(data), sizeData / 2400)

    params, cov = curve_fit(bimodal, edges, heights, expected)

    plt.plot(edges, bimodal(edges, *params), color='red', lw=3, label='model')
    plt.plot(edges, heights)
    plt.show()


def folowMaxValue(vol, start_seed, direction_seed, stop_distance):
    diameter = []
    distance = []
    c_s = [start_seed[2], start_seed[1], start_seed[0]]
    d_s = [direction_seed[2], direction_seed[1], direction_seed[0]]

    vol2 = np.zeros((vol.shape[0], vol.shape[1], vol.shape[2]))

    last_seeds = []
    last_seeds.append(c_s)
    last_seeds.append(d_s)

    diameter.append(vol[c_s[0], c_s[1], c_s[2]])
    distance.append(0)
    vol[c_s[0], c_s[1], c_s[2]] = 0
    vol2[c_s[0], c_s[1], c_s[2]] = 1.0

    old_dir = np.zeros((10, 3))
    old_dir[0, :] = [d_s[0] - c_s[0], d_s[1] - c_s[1], d_s[2] - c_s[2]]

    dicToReturn = {'Image': [], 'Distance': [], 'Diameter': []}

    for index in range(1, stop_distance):
        x1 = c_s[0] - 1
        x2 = c_s[0] + 2
        y1 = c_s[1] - 1
        y2 = c_s[1] + 2
        z1 = c_s[2] - 1
        z2 = c_s[2] + 2

        if x1 == (-1):
            x1 = 0
        if x2 == (vol.shape[0] + 1):
            x2 = vol.shape[0]

        if y1 == (-1):
            y1 = 0
        if y2 == (vol.shape[1] + 1):
            y2 = vol.shape[1]

        if z1 == (-1):
            z1 = 0
        if z2 == (vol.shape[2] + 1):
            z2 = vol

    for z in range(0, vol.shape[0]):
        looking = vol[x1:x2, y1:y2, z1:z2]

        size = 1000
        while (size > 4):
            mean_value = np.mean(np.mean(looking[np.nonzero(looking)]))
            looking[looking < mean_value] = 0
            coord = np.transpose(np.nonzero(looking))
            size = len(coord)

        if not coord.size:
            break
        else:
            error2Beat = 100000
            distance.append(0.0)
            diameter.append(0.0)

            for c in coord:
                mean_dir = [np.mean(old_dir[:, 0]), np.mean(old_dir[:, 1]), np.mean(old_dir[:, 2])]
                c_dir = [c[0] + x1 - c_s[0], c[1] + y1 - c_s[1], c[2] + z1 - c_s[2]]
                error_dir = np.sqrt(
                    (c_dir[0] - mean_dir[0]) ** 2 + (c_dir[1] - mean_dir[1]) ** 2 + (c_dir[2] - mean_dir[2]) ** 2)
                error_value = np.sqrt((vol[c[0] + x1, c[1] + y1, c[2] + z1] - vol[c_s[0], c_s[1], c_s[2]]) ** 2) / 10.0
                error = error_dir + error_value
                if error < error2Beat:
                    c_choose = c
                    distance[-1] = distance[-2] + np.sqrt(
                        (c_s[0] - (c_choose[0] + x1)) ** 2 + (c_s[1] - (c_choose[1] + y1)) ** 2 + (
                                    c_s[2] - (c_choose[2] + z1)) ** 2)
                    diameter[-1] = vol[c_choose[0] + x1, c_choose[1] + y1, c_choose[2] + z1]

                    vol[c[0] + x1, c[1] + y1, c[2] + z1] = 0

                old_dir[index % 10] = [c_choose[0] + x1 - c_s[0], c_choose[1] + y1 - c_s[1], c_choose[2] + z1 - c_s[2]]

                c_s = [c_choose[0] + x1, c_choose[1] + y1, c_choose[2] + z1]
                vol[c_s[0], c_s[1], c_s[2]] = 0
                vol2[c_s[0], c_s[1], c_s[2]] = index
                index += 1

    dicToReturn['Image'] = vol2
    dicToReturn['Diameter'] = diameter
    dicToReturn['Distance'] = distance

    return dicToReturn


def folowLine(vol, start_seed, direction_seed, stop_seed):
    diameter = []
    distance = []
    c_s = [start_seed[2], start_seed[1], start_seed[0]]
    d_s = [direction_seed[2], direction_seed[1], direction_seed[0]]

    vol2 = np.zeros((vol.shape[0], vol.shape[1], vol.shape[2]))

    last_seeds = []
    last_seeds.append(c_s)
    last_seeds.append(d_s)

    diameter.append(vol[c_s[0], c_s[1], c_s[2]])
    distance.append(0)
    vol[c_s[0], c_s[1], c_s[2]] = 0

    old_dir = np.zeros((5, 3))
    old_dir[0, :] = [d_s[0] - c_s[0], d_s[1] - c_s[1], d_s[2] - c_s[2]]

    index = 1
    dicToReturn = {'Image': [], 'Distance': [], 'Diameter': []}
    while (not np.array_equal(c_s, [stop_seed[2], stop_seed[1], stop_seed[0]])):

        vol2[c_s[0], c_s[1], c_s[2]] = 1.0

        x1 = c_s[0] - 1
        x2 = c_s[0] + 2
        y1 = c_s[1] - 1
        y2 = c_s[1] + 2
        z1 = c_s[2] - 1
        z2 = c_s[2] + 2

        if x1 == (-1):
            x1 = 0
        if x2 == (vol.shape[0] + 1):
            x2 = vol.shape[0]

        if y1 == (-1):
            y1 = 0
        if y2 == (vol.shape[1] + 1):
            y2 = vol.shape[1]

        if z1 == (-1):
            z1 = 0
        if z2 == (vol.shape[2] + 1):
            z2 = vol.shape[2]

        looking = vol[x1:x2, y1:y2, z1:z2]
        coord = np.transpose(np.nonzero(looking))
        if not coord.size:
            break
        else:
            error2Beat = 100000
            for c in coord:
                mean_dir = [np.mean(old_dir[:, 0]), np.mean(old_dir[:, 1]), np.mean(old_dir[:, 2])]
                c_dir = [c[0] + x1 - c_s[0], c[1] + y1 - c_s[1], c[2] + z1 - c_s[2]]
                error = (c_dir[0] - mean_dir[0]) ** 2 + (c_dir[1] - mean_dir[1]) ** 2 + (c_dir[2] - mean_dir[2]) ** 2
                if error < error2Beat:
                    c_choose = c

            old_dir[index % 5] = [c_choose[0] + x1 - c_s[0], c_choose[1] + y1 - c_s[1], c_choose[2] + z1 - c_s[2]]
            distance.append(distance[-1] + np.sqrt(
                (c_s[0] - (c_choose[0] + x1)) ** 2 + (c_s[1] - (c_choose[1] + y1)) ** 2 + (
                            c_s[2] - (c_choose[2] + z1)) ** 2))
            diameter.append(vol[c_choose[0] + x1, c_choose[1] + y1, c_choose[2] + z1])
            c_s = [c_choose[0] + x1, c_choose[1] + y1, c_choose[2] + z1]
            vol[c_s[0], c_s[1], c_s[2]] = 0
            index += 1

    dicToReturn['Image'] = vol2
    dicToReturn['Diameter'] = diameter
    dicToReturn['Distance'] = distance

    return dicToReturn


def vtk_import_numpy_array(np_array):
    np_array = (255 * (np_array - np.min(np_array))) / (np.max(np_array) - np.min(np_array))
    np_array = np_array.astype(np.uint8)
    np_array[np_array < 0] = 0
    np_array[np_array > 255] = 255

    data_importer = vtk.vtkImageImport()
    data_importer.CopyImportVoidPointer(np_array, np_array.nbytes)
    data_importer.SetDataScalarTypeToUnsignedChar()
    data_importer.SetNumberOfScalarComponents(1)
    data_importer.SetDataExtent(0, np_array.shape[2] - 1, 0, np_array.shape[1] - 1, 0, np_array.shape[0] - 1)
    data_importer.SetWholeExtent(0, np_array.shape[2] - 1, 0, np_array.shape[1] - 1, 0, np_array.shape[0] - 1)

    return data_importer


def importSTL(filename):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()

    outputVol = np.zeros((500, 500, 500))
    outputVTK = vtk_import_numpy_array(outputVol)

    dataToStencil = vtk.vtkPolyDataToImageStencil()
    dataToStencil.SetInputConnection(reader.GetOutputPort())

    dataToStencil.SetOutputSpacing(0.5, 0.5, 0.5)
    dataToStencil.SetOutputOrigin(0.0, 0.0, 0.0)
    dataToStencil.Update()

    stencil = vtk.vtkImageStencil()
    stencil.SetInput(outputVTK.GetOutput())
    stencil.SetStencil(dataToStencil.GetOutput())
    stencil.ReverseStencilOn()
    stencil.Update()

    im = stencil.GetOutput()
    rows, cols, lines = im.GetDimensions()
    sc = im.GetPointData().GetScalars()

    a = vtk_to_numpy(sc)
    a = a.reshape((lines, cols, rows))
    return a


def smooth_3D(vol, nbIter, Factor):
    vtk_importer = vtk_import_numpy_array(vol)

    threshold = vtk.vtkImageThreshold()
    threshold.SetInputConnection(vtk_importer.GetOutputPort())
    threshold.ThresholdByLower(125)
    threshold.ReplaceInOn()
    threshold.SetInValue(0)
    threshold.ReplaceOutOn()
    threshold.SetOutValue(1)
    threshold.Update()

    dmc = vtk.vtkDiscreteMarchingCubes()
    dmc.SetInputConnection(threshold.GetOutputPort())
    dmc.GenerateValues(1, 1, 1)
    dmc.Update()

    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputConnection(dmc.GetOutputPort())
    smoother.SetNumberOfIterations(nbIter)
    smoother.BoundarySmoothingOff()

    smoother.NonManifoldSmoothingOn()
    smoother.NormalizeCoordinatesOn()
    smoother.GenerateErrorScalarsOn()

    smoother.Update()
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName('./Data/' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + '.stl')
    stlWriter.SetInputConnection(smoother.GetOutputPort())
    stlWriter.Write()

    dmc = vtk.vtkPolyDataNormals()
    dmc.SetInputConnection(smoother.GetOutputPort())

    outputVol = np.zeros((vol.shape[0] * Factor, vol.shape[1] * Factor, vol.shape[2] * Factor))
    outputVTK = vtk_import_numpy_array(outputVol)

    dataToStencil = vtk.vtkPolyDataToImageStencil()
    dataToStencil.SetInputConnection(dmc.GetOutputPort())

    dataToStencil.SetOutputSpacing(1.0 / Factor, 1.0 / Factor, 1.0 / Factor)
    dataToStencil.SetOutputOrigin(0.0, 0.0, 0.0)
    dataToStencil.Update()

    stencil = vtk.vtkImageStencil()
    stencil.SetInput(outputVTK.GetOutput())
    stencil.SetStencil(dataToStencil.GetOutput())
    stencil.ReverseStencilOn()
    stencil.Update()

    im = stencil.GetOutput()
    rows, cols, lines = im.GetDimensions()
    sc = im.GetPointData().GetScalars()

    a = vtk_to_numpy(sc)
    a = a.reshape((lines, cols, rows))

    return a


def dwt_denoise(vol, waveletType, levels, alpha):
    wavelet = pywt.Wavelet(waveletType)

    X = vol.shape[1] - 1
    Y = vol.shape[2] - 1

    NewImage = np.zeros((vol.shape[0], vol.shape[1], vol.shape[2]))
    levels = int(np.floor(np.log2(vol.shape[1])))
    th = alpha * np.sqrt(2 * np.log2(vol.shape[1]))

    for z in range(0, vol.shape[0]):
        image = vol[z, :, :]
        WaveletCoeffs = pywt.wavedec2(image, wavelet, level=levels)
        NewWaveletCoeffs = map(lambda x: pywt.thresholding.soft(x, th), WaveletCoeffs)
        Image = pywt.waverec2(NewWaveletCoeffs, wavelet)
        NewImage[z, :X, :Y] = Image[:X, :Y]

    return NewImage


""" Image Processing """


def mutualInformation(ImageStack, TimeRef, FlagROI, ROI):
    MutualInfo = []

    for z in range(TimeRef[0] + 1, TimeRef[1]):

        if FlagROI:
            image1 = ImageStack[z - 1, ROI[0]:ROI[1], ROI[2]:ROI[3]]
            image2 = ImageStack[z, ROI[0]:ROI[1], ROI[2]:ROI[3]]
        else:
            image1 = ImageStack[z - 1, :, :]
            image2 = ImageStack[z, :, :]

        c = ee.vectorize(image1.flatten())
        d = ee.vectorize(image2.flatten())

        mi = ee.mi(c, d)
        MutualInfo.append(mi)

    return MutualInfo


def giveChessImage(Im1, Im2, BlocSpeed):
    SizeZVol = Im1.shape[0]
    SizeXVol = Im1.shape[1]
    SizeYVol = Im1.shape[2]

    iterx = 0
    itery = 0
    iterz = 0

    chessMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))

    for xRef in np.arange(0, SizeXVol - 1, BlocSpeed):
        itery = 0
        iterx += 1
        for yRef in np.arange(0, SizeYVol - 1, BlocSpeed):
            iterz = 0
            itery += 1
            for zRef in np.arange(0, SizeZVol - 1, BlocSpeed):
                iterz += 1

                xMinRef = xRef - int(BlocSpeed / 2)
                yMinRef = yRef - int(BlocSpeed / 2)
                zMinRef = zRef - int(BlocSpeed / 2)

                xMaxRef = xRef + int(BlocSpeed / 2)
                yMaxRef = yRef + int(BlocSpeed / 2)
                zMaxRef = zRef + int(BlocSpeed / 2)

                if xMinRef < 0:
                    xMinRef = 0
                if yMinRef < 0:
                    yMinRef = 0
                if zMinRef < 0:
                    zMinRef = 0

                if xMaxRef >= SizeXVol:
                    xMaxRef = SizeXVol - 1
                if yMaxRef >= SizeYVol:
                    yMaxRef = SizeYVol - 1
                if zMaxRef >= SizeZVol:
                    zMaxRef = SizeZVol - 1

                if ((iterz + itery + iterx) % 2) == 0:
                    chessMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] = Im1[zMinRef:zMaxRef, xMinRef:xMaxRef,
                                                                                  yMinRef:yMaxRef]
                else:
                    chessMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] = Im2[zMinRef:zMaxRef, xMinRef:xMaxRef,
                                                                                  yMinRef:yMaxRef]

    return chessMap


def func(x, a, b):
    return a * (1.0 - np.exp(-x / b))


def resizeImage(image2resample, interpolator, spacing):
    img = imageFromNumpyToITK(image2resample)

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

    image = imageFromITKToNumpy(img)

    return image


def ComputeMaskLocalDensity(mask, ws, ov):
    SizeZVol = mask.shape[0]
    SizeXVol = mask.shape[1]
    SizeYVol = mask.shape[2]

    BlocSpeed = int((1 - ov) * ws)
    if BlocSpeed == 0:
        print("The overlap ration is to big the bloc matching speed is set to 1 pixel")
        BlocSpeed = 1

    densityMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))
    iterationMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))

    if SizeZVol == 1:
        BlocSpeedZ = 1
    else:
        BlocSpeedZ = BlocSpeed

    for xRef in np.arange(0, SizeXVol, BlocSpeed):
        for yRef in np.arange(0, SizeYVol, BlocSpeed):
            for zRef in np.arange(0, SizeZVol, BlocSpeedZ):

                xMinRef = xRef - int(ws / 2)
                yMinRef = yRef - int(ws / 2)
                zMinRef = zRef - int(ws / 2)

                xMaxRef = xRef + int(ws / 2)
                yMaxRef = yRef + int(ws / 2)
                zMaxRef = zRef + int(ws / 2)

                if xMinRef < 0:
                    xMinRef = 0
                if yMinRef < 0:
                    yMinRef = 0
                if zMinRef < 0:
                    zMinRef = 0

                if xMaxRef >= SizeXVol:
                    xMaxRef = SizeXVol
                if yMaxRef >= SizeYVol:
                    yMaxRef = SizeYVol
                if zMaxRef >= SizeZVol:
                    zMaxRef = SizeZVol

                BlocRef = mask[int((zMinRef + zMaxRef) / 2), xMinRef:xMaxRef, yMinRef:yMaxRef]
                #
                iterationMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] += 1
                all_labels = measure.label(BlocRef, background=0)
                all_labels[all_labels == -1] = 0

                densityMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] += float(np.max(all_labels, axis=None))

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(densityMap, iterationMap)
        c[c == np.inf] = 0
        c = np.nan_to_num(c)

    return c


def ComputeMaskLocalVolumeDensity(mask, ws, ov):
    SizeZVol = mask.shape[0]
    SizeXVol = mask.shape[1]
    SizeYVol = mask.shape[2]

    BlocSpeed = int((1 - ov) * ws)
    if BlocSpeed == 0:
        print("The overlap ration is to big the bloc matching speed is set to 1 pixel")
        BlocSpeed = 1

    densityMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))
    iterationMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))

    if SizeZVol == 1:
        BlocSpeedZ = 1
    else:
        BlocSpeedZ = BlocSpeed

    for xRef in np.arange(0, SizeXVol, BlocSpeed):
        for yRef in np.arange(0, SizeYVol, BlocSpeed):
            for zRef in np.arange(0, SizeZVol, BlocSpeedZ):

                xMinRef = xRef - int(ws / 2)
                yMinRef = yRef - int(ws / 2)
                zMinRef = zRef - int(ws / 2)

                xMaxRef = xRef + int(ws / 2)
                yMaxRef = yRef + int(ws / 2)
                zMaxRef = zRef + int(ws / 2)

                if xMinRef < 0:
                    xMinRef = 0
                if yMinRef < 0:
                    yMinRef = 0
                if zMinRef < 0:
                    zMinRef = 0

                if xMaxRef >= SizeXVol:
                    xMaxRef = SizeXVol
                if yMaxRef >= SizeYVol:
                    yMaxRef = SizeYVol
                if zMaxRef >= SizeZVol:
                    zMaxRef = SizeZVol

                BlocRef = mask[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef]

                densityMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] += (
                            np.sum(BlocRef, axis=None) / (BlocRef.shape[0] * BlocRef.shape[1] * BlocRef.shape[2]))
                iterationMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(densityMap, iterationMap)
        c[c == np.inf] = 0
        c = np.nan_to_num(c)

    return c


def ComputeMaskLocalNumberDensity(mask, ws, ov):
    pass


def ComputeMaskLocalSurfaceDensity(mask, ws, ov):
    pass


def ComputeDensityChange(hpImage, lpImage, ws, ov):
    SizeZVol = hpImage.shape[0]
    SizeXVol = hpImage.shape[1]
    SizeYVol = hpImage.shape[2]

    hpImage = (254.0) / (1.0 + np.exp(-(hpImage - 200.0) / -100.0)) + 2.0
    lpImage = (254.0) / (1.0 + np.exp(-(lpImage - 200.0) / -100.0)) + 2.0

    BlocSpeed = int((1 - ov) * ws)
    if BlocSpeed == 0:
        print("The overlap ration is to big the bloc matching speed is set to 1 pixel")
        BlocSpeed = 1

    densityMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))
    iterationMap = np.zeros((SizeZVol, SizeXVol, SizeYVol))

    'Ready For Loop'

    for xRef in np.arange(0, SizeXVol - 1, BlocSpeed):
        for yRef in np.arange(0, SizeYVol - 1, BlocSpeed):
            for zRef in np.arange(0, SizeZVol - 1, BlocSpeed):

                xMinRef = xRef - int(ws / 2)
                yMinRef = yRef - int(ws / 2)
                zMinRef = zRef - int(ws / 2)

                xMaxRef = xRef + int(ws / 2)
                yMaxRef = yRef + int(ws / 2)
                zMaxRef = zRef + int(ws / 2)

                if xMinRef < 0:
                    xMinRef = 0
                if yMinRef < 0:
                    yMinRef = 0
                if zMinRef < 0:
                    zMinRef = 0

                if xMaxRef >= SizeXVol:
                    xMaxRef = SizeXVol - 1
                if yMaxRef >= SizeYVol:
                    yMaxRef = SizeYVol - 1
                if zMaxRef >= SizeZVol:
                    zMaxRef = SizeZVol - 1

                BlocRefHP = hpImage[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] + 1.0
                BlocRefLP = lpImage[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] + 1.0
                BlocRef = (BlocRefLP - BlocRefHP)

                densityMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] += stats.nanmedian(BlocRef, axis=None)
                iterationMap[zMinRef:zMaxRef, xMinRef:xMaxRef, yMinRef:yMaxRef] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(densityMap, iterationMap)
        c[c == np.inf] = 0
        c = np.nan_to_num(c)

    return c


def computeVentilationMap(Im1, OverlapRatio, SizeBloc, startingImage, time_step):
    SizeZVol = Im1.shape[0]
    SizeXVol = Im1.shape[1]
    SizeYVol = Im1.shape[2]

    BlocSpeed = int((1 - OverlapRatio) * SizeBloc)
    if BlocSpeed == 0:
        print("The overlap ration is to big the bloc matching speed is set to 1 pixel")
        BlocSpeed = 1

    ventilationMap = np.zeros((1, SizeXVol, SizeYVol))
    IterationMap = np.zeros((1, SizeXVol, SizeYVol))
    t_array = np.arange(0, SizeZVol * time_step, time_step)

    for xRef in np.arange(0, SizeXVol - 1, BlocSpeed):
        for yRef in np.arange(0, SizeYVol - 1, BlocSpeed):
            xMinRef = xRef - int(SizeBloc / 2)
            yMinRef = yRef - int(SizeBloc / 2)

            xMaxRef = xRef + int(SizeBloc / 2)
            yMaxRef = yRef + int(SizeBloc / 2)

            if xMinRef < 0:
                xMinRef = 0
            if yMinRef < 0:
                yMinRef = 0

            if xMaxRef >= SizeXVol:
                xMaxRef = SizeXVol - 1
            if yMaxRef >= SizeYVol:
                yMaxRef = SizeYVol - 1

            BlocRef = Im1[:, xMinRef:xMaxRef, yMinRef:yMaxRef]

            c_array = np.array(np.median(np.median(BlocRef, axis=1), axis=1))
            if np.median(c_array) != 0:
                try:
                    popt, pcov = curve_fit(func, t_array, c_array)
                    if not np.isnan(popt[1]):
                        ventilationMap[:, xMinRef:xMaxRef, yMinRef:yMaxRef] += popt[1]
                        IterationMap[:, xMinRef:xMaxRef, yMinRef:yMaxRef] += 1

                except:
                    pass

    ventilationMap = np.nan_to_num(ventilationMap / IterationMap)
    ventilationMap[ventilationMap > 5.0] = 5.0
    ventilationMap[ventilationMap < 0] = 0
    return ventilationMap


def equalize(image, nzone, zones):
    toReturn = np.zeros((image.shape[0], image.shape[1], image.shape[2]))

    x11 = zones[0][0]
    y11 = zones[0][1]
    z11 = zones[0][2]
    x21 = zones[0][3]
    y21 = zones[0][4]
    z21 = zones[0][5]

    if x11 == x21:
        dir_i = 0
        meanZone1 = np.mean(image[x11, y11:y21, z11:z21])
    if y11 == y21:
        dir_i = 1
        meanZone1 = np.mean(image[x11:x21, y11, z11:z21])
    if z11 == z21:
        dir_i = 2
        meanZone1 = np.mean(image[x11:x21, y11:y21, z11])

    if (nzone == 2):
        x12 = zones[1][0]
        y12 = zones[1][1]
        z12 = zones[1][2]
        x22 = zones[1][3]
        y22 = zones[1][4]
        z22 = zones[1][5]

        if (x12 == x22) and (dir_i == 0):
            meanZone2 = np.mean(image[x12, y12:y22, z12:z22])
        elif y12 == y22 and (dir_i == 1):
            meanZone2 = np.mean(image[x12:x22, y12, z12:z22])
        elif z12 == z22 and (dir_i == 2):
            meanZone2 = np.mean(image[x12:x22, y12:y22, z12])

    for i in range(image.shape[dir_i]):

        if dir_i == 0:
            currentSlice = image[i, :, :]
            mean_zone1_i = np.mean(currentSlice[y11:y21, z11:z21])
            if nzone == 1:
                toReturn[i, :, :] = currentSlice + meanZone1 - mean_zone1_i
            elif nzone == 2:
                mean_zone2_i = np.mean(currentSlice[y12:y22, z12:z22])
                toReturn[i, :, :] = (currentSlice - mean_zone2_i) * (meanZone1 - meanZone2) / (
                            mean_zone1_i - mean_zone2_i) + meanZone2

        elif dir_i == 1:
            currentSlice = image[:, i, :]
            mean_zone1_i = np.mean(currentSlice[x11:x21, z11:z21])
            if nzone == 1:
                toReturn[:, i, :] = currentSlice + meanZone1 - mean_zone1_i
            elif nzone == 2:
                mean_zone2_i = np.mean(currentSlice[x12:x22, z12:z22])
                toReturn[:, i, :] = (currentSlice - mean_zone2_i) * (meanZone1 - meanZone2) / (
                            mean_zone1_i - mean_zone2_i) + meanZone2

        elif dir_i == 2:
            currentSlice = image[:, :, i]
            mean_zone1_i = np.mean(currentSlice[x11:x21, y11:y21])
            if nzone == 1:
                toReturn[:, :, i] = currentSlice + meanZone1 - mean_zone1_i
            elif nzone == 2:
                mean_zone2_i = np.mean(currentSlice[x12:x22, y12:y22])
                toReturn[:, :, i] = (currentSlice - mean_zone2_i) * (meanZone1 - meanZone2) / (
                            mean_zone1_i - mean_zone2_i) + meanZone2

    return toReturn


def Threshold(vol, min_th, max_th, inside_value=1, outside_value=0):
    ITK_Vol = imageFromNumpyToITK(vol)
    binaryThresh = sitk.BinaryThresholdImageFilter()

    binaryThresh.SetLowerThreshold(min_th)
    binaryThresh.SetUpperThreshold(max_th)
    binaryThresh.SetInsideValue(inside_value)
    binaryThresh.SetOutsideValue(outside_value)
    ITK_Vol = binaryThresh.Execute(ITK_Vol)

    return imageFromITKToNumpy(ITK_Vol)


def SegWatershed(inputDataToSeg, level, markLine, connected):
    ITK_Vol = imageFromNumpyToITK(inputDataToSeg)
    watershedF = sitk.MorphologicalWatershedImageFilter()
    watershedF.SetLevel(level)
    watershedF.SetFullyConnected(connected)
    watershedF.SetMarkWatershedLine(markLine)

    ITK_Vol = watershedF.Execute(ITK_Vol)
    image = imageFromITKToNumpy(ITK_Vol)

    return image


def SegConnectedThreshold(vol, val_min, val_max, seedListToSegment):
    ITK_Vol = imageFromNumpyToITK(vol)

    segmentationFilter = sitk.ConnectedThresholdImageFilter()

    for seed in seedListToSegment:
        seedItk = (seed[0], seed[1], seed[2])
        segmentationFilter.AddSeed(seedItk)

    segmentationFilter.SetLower(val_min)
    segmentationFilter.SetUpper(val_max)
    segmentationFilter.SetReplaceValue(1)
    ITK_Vol = segmentationFilter.Execute(ITK_Vol)
    image = imageFromITKToNumpy(ITK_Vol)

    return image.astype(np.uint8)


def rotation(image, rotation):
    shapeNewImage = ndimage.rotate(image[0, :, :], rotation, reshape=True).shape
    newImage = np.zeros((image.shape[0], shapeNewImage[0], shapeNewImage[1]))
    for z in range(0, image.shape[0]):
        newImage[z, :, :] = ndimage.rotate(image[z, :, :], rotation, reshape=True)

    return newImage


def histogram(array, minBin, maxBin, nbBin):
    bins = np.arange(minBin, maxBin, (maxBin - minBin) / nbBin)
    histogram = np.histogram(array, bins)

    return histogram


def SegConnectedThresholdC(vol, radius, multi, iterN, seedListToSegment):
    ITK_Vol = imageFromNumpyToITK(vol)

    segmentationFilter = sitk.ConfidenceConnectedImageFilter()

    for seed in seedListToSegment:
        seedItk = (seed[0], seed[1], seed[2])
        segmentationFilter.AddSeed(seedItk)

    segmentationFilter.SetInitialNeighborhoodRadius(radius)
    segmentationFilter.SetMultiplier(multi)
    segmentationFilter.SetNumberOfIterations(iterN)
    segmentationFilter.SetReplaceValue(1)
    ITK_Vol = segmentationFilter.Execute(ITK_Vol)
    image = imageFromITKToNumpy(ITK_Vol)
    return image.astype(np.uint8)


def FastMarching(vol, stopValue, seedListToSegment):
    ITK_Vol = imageFromNumpyToITK(vol)

    fast_marching = sitk.FastMarchingImageFilter()

    for seed in seedListToSegment:
        seedItk = (seed[0], seed[1], seed[2])
        fast_marching.AddTrialPoint(seedItk)

    fast_marching.SetStoppingValue(stopValue)
    ITK_Vol = fast_marching.Execute(ITK_Vol)

    return imageFromITKToNumpy(ITK_Vol)


def ShapeDetectionLS(vol_LS, vol_EP, maxRMSError, propaScaling, curvScaling, nbIter, radius, seedListToSegment):
    vol_EP = vol_EP.astype(float)
    ITK_EP = imageFromNumpyToITK(vol_EP)

    if radius != 0.0:

        vol_LS = np.zeros((vol_EP.shape[0], vol_EP.shape[1], vol_EP.shape[2]))

        for seed in seedListToSegment:

            start_i, end_i = seed[0] - int(radius) - 1, seed[0] + int(radius) + 1
            start_j, end_j = seed[1] - int(radius) - 1, seed[1] + int(radius) + 1
            start_k, end_k = seed[2] - int(radius) - 1, seed[2] + int(radius) + 1

            for i in range(start_i, end_i):
                for j in range(start_j, end_j):
                    for k in range(start_k, end_k):
                        distance = int(((i - seed[0]) ** 2 + (j - seed[1]) ** 2 + (k - seed[2]) ** 2) ** 0.5)
                        if distance <= radius:
                            vol_LS[k, j, i] = 1.0

    if vol_LS.max() > 1.0:
        vol_LS /= vol_LS.max()

    vol_LS *= -1.0
    vol_LS += 0.5

    ITK_LS = imageFromNumpyToITK(vol_LS)
    shapeDetect = sitk.ShapeDetectionLevelSetImageFilter()
    shapeDetect.SetMaximumRMSError(maxRMSError)
    shapeDetect.SetPropagationScaling(propaScaling)
    shapeDetect.SetCurvatureScaling(curvScaling)
    shapeDetect.SetNumberOfIterations(nbIter)
    ITK_Vol = shapeDetect.Execute(ITK_LS, ITK_EP)

    return imageFromITKToNumpy(ITK_Vol)


def GeoDetectionLS(vol_LS, vol_EP, maxRMSError, propaScaling, curvScaling, advScaling, nbIter, radius,
                   seedListToSegment):
    vol_EP = vol_EP.astype(float)
    ITK_EP = imageFromNumpyToITK(vol_EP)

    if radius != 0.0:
        vol_LS = np.zeros((vol_EP.shape[0], vol_EP.shape[1], vol_EP.shape[2]))

        for seed in seedListToSegment:

            start_i, end_i = seed[0] - int(radius) - 1, seed[0] + int(radius) + 1
            start_j, end_j = seed[1] - int(radius) - 1, seed[1] + int(radius) + 1
            start_k, end_k = seed[2] - int(radius) - 1, seed[2] + int(radius) + 1

            for i in range(start_i, end_i):
                for j in range(start_j, end_j):
                    for k in range(start_k, end_k):
                        distance = int(((i - seed[0]) ** 2 + (j - seed[1]) ** 2 + (k - seed[2]) ** 2) ** 0.5)
                        if distance <= radius:
                            vol_LS[k, j, i] = 1.0

    if vol_LS.max() > 1.0:
        vol_LS /= vol_LS.max()

    vol_LS *= -1.0
    vol_LS += 0.5

    ITK_LS = imageFromNumpyToITK(vol_LS)
    shapeDetect = sitk.GeodesicActiveContourLevelSetImageFilter()
    shapeDetect.SetMaximumRMSError(maxRMSError)
    shapeDetect.SetPropagationScaling(propaScaling)
    shapeDetect.SetCurvatureScaling(curvScaling)
    shapeDetect.SetAdvectionScaling(advScaling)
    shapeDetect.SetNumberOfIterations(nbIter)
    ITK_Vol = shapeDetect.Execute(ITK_LS, ITK_EP)

    return imageFromITKToNumpy(ITK_Vol)


def InterpolateDataPoints(dataPoints, ImageSize):
    ImageD0 = np.zeros((ImageSize[0], ImageSize[1], ImageSize[2]))
    contoursD0 = dataPoints['Direction0']

    ImageD1 = np.zeros((ImageSize[0], ImageSize[1], ImageSize[2]))
    contoursD1 = dataPoints['Direction1']

    ImageD2 = np.zeros((ImageSize[0], ImageSize[1], ImageSize[2]))
    contoursD2 = dataPoints['Direction2']

    ImageDT = np.zeros((ImageSize[0], ImageSize[1], ImageSize[2]))
    all_points = []

    if any(contoursD0):
        new_arr = []
        for i, curve in enumerate(contoursD0):
            if any(curve):
                if i == 0:
                    polygon = []
                    for point in curve:
                        polygon.append((point[1], point[0]))
                    img = Image.new('L', (ImageSize[1], ImageSize[2]), 0)
                    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                    mask = np.array(img)

                    ImageD0[:curve[0][2]] = mask.T
                elif i == (len(contoursD0) - 1):
                    polygon = []
                    for point in curve:
                        polygon.append((point[1], point[0]))

                    img = Image.new('L', (ImageSize[1], ImageSize[2]), 0)
                    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                    mask = np.array(img)
                    ImageD0[curve[0][2]:] = mask.T

                for coord in curve:
                    if len(coord) != 0:
                        new_arr.append(coord)
                        all_points.append(coord)

        nx = ImageSize[1]
        ny = ImageSize[2]
        nz = ImageSize[0]

        tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
        u_new = np.linspace(u.min(), u.max(), len(coord) * 100)
        x, y, z = splev(u_new, tck, der=0)

        new_arr = np.array(new_arr)

        g = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))
        positions = np.vstack(map(np.ravel, g))

        hull = scipy.spatial.Delaunay(new_arr)

        mask = (hull.find_simplex(positions.T) >= 0)

        positions = positions[:, mask]
        positions = positions[::-1]

        ImageD0[positions[0, :], positions[1, :], positions[2, :]] = 1
    else:
        ImageD0[:, :, :] = 1

    if any(contoursD1):
        new_arr = []
        for i, curve in enumerate(contoursD1):
            if any(curve):
                if i == 0:
                    polygon = []
                    for point in curve:
                        polygon.append((point[2], point[0]))

                    img = Image.new('L', (ImageSize[0], ImageSize[2]), 0)
                    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                    mask = np.array(img)

                    for x_i in range(0, int(curve[0][1])):
                        ImageD1[:, x_i, :] = mask.T

                elif i == (len(contoursD1) - 1):
                    polygon = []
                    for point in curve:
                        polygon.append((point[2], point[0]))

                    img = Image.new('L', (ImageSize[0], ImageSize[2]), 0)
                    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                    mask = np.array(img)
                    for x_i in range(int(curve[0][1]), ImageSize[1]):
                        ImageD1[:, x_i, :] = mask.T

                for coord in curve:
                    if len(coord) != 0:
                        new_arr.append(coord)
                        all_points.append(coord)

        nx = ImageSize[1]
        ny = ImageSize[2]
        nz = ImageSize[0]

        tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
        u_new = np.linspace(u.min(), u.max(), len(coord) * 300)
        x, y, z = splev(u_new, tck, der=0)

        new_arr = np.array(new_arr)

        g = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))
        positions = np.vstack(map(np.ravel, g))

        hull = scipy.spatial.Delaunay(new_arr)

        mask = (hull.find_simplex(positions.T) >= 0)

        positions = positions[:, mask]
        positions = positions[::-1]

        ImageD1[positions[0, :], positions[1, :], positions[2, :]] = 1
    else:
        ImageD1[:, :, :] = 1

    if any(contoursD2):
        new_arr = []
        for i, curve in enumerate(contoursD2):
            if any(curve):
                if i == 0:
                    polygon = []
                    for point in curve:
                        polygon.append((point[2], point[1]))

                    img = Image.new('L', (ImageSize[0], ImageSize[1]), 0)
                    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                    mask = np.array(img)

                    for y_i in range(0, int(curve[0][0])):
                        ImageD2[:, :, y_i] = mask.T

                elif i == (len(contoursD2) - 1):
                    polygon = []
                    for point in curve:
                        polygon.append((point[2], point[1]))

                    img = Image.new('L', (ImageSize[0], ImageSize[1]), 0)
                    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                    mask = np.array(img)
                    for y_i in range(int(curve[0][0]), ImageSize[2]):
                        ImageD2[:, :, y_i] = mask.T

                for coord in curve:
                    if len(coord) != 0:
                        new_arr.append(coord)
                        all_points.append(coord)

        nx = ImageSize[1]
        ny = ImageSize[2]
        nz = ImageSize[0]

        tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
        u_new = np.linspace(u.min(), u.max(), len(coord) * 300)
        x, y, z = splev(u_new, tck, der=0)

        new_arr = np.array(new_arr)

        g = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))
        positions = np.vstack(map(np.ravel, g))

        hull = scipy.spatial.Delaunay(new_arr)

        mask = (hull.find_simplex(positions.T) >= 0)

        positions = positions[:, mask]
        positions = positions[::-1]

        ImageD2[positions[0, :], positions[1, :], positions[2, :]] = 1
    else:
        ImageD2[:, :, :] = 1

    nx = ImageSize[1]
    ny = ImageSize[2]
    nz = ImageSize[0]

    tck, u = splprep(np.array(new_arr).T, u=None, s=0.0, per=1)
    u_new = np.linspace(u.min(), u.max(), len(coord) * 300)
    x, y, z = splev(u_new, tck, der=0)

    new_arr = np.array(all_points)

    g = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))
    positions = np.vstack(map(np.ravel, g))

    hull = scipy.spatial.Delaunay(new_arr)

    mask = (hull.find_simplex(positions.T) >= 0)

    positions = positions[:, mask]
    positions = positions[::-1]

    ImageDT[positions[0, :], positions[1, :], positions[2, :]] = 1

    return np.multiply(np.multiply(np.multiply(ImageD2, ImageD1), ImageD0), ImageDT)


def anisotropic_diffusion(vol, time_step, conductance, nbIter):
    ITK_Vol = imageFromNumpyToITK(vol)

    curvDiff = sitk.CurvatureAnisotropicDiffusionImageFilter()
    curvDiff.SetTimeStep(time_step)
    curvDiff.SetConductanceParameter(conductance)
    curvDiff.SetNumberOfIterations(nbIter)
    ITK_Vol = curvDiff.Execute(ITK_Vol)

    return imageFromITKToNumpy(ITK_Vol)


def zero_crossing(vol, var, maxErr):
    ITK_Vol = imageFromNumpyToITK(vol)
    zeroCrossingEdge = sitk.ZeroCrossingBasedEdgeDetectionImageFilter()
    zeroCrossingEdge.SetVariance(var)
    zeroCrossingEdge.SetMaximumError(maxErr)
    ITK_Vol = zeroCrossingEdge.Execute(ITK_Vol)
    return imageFromITKToNumpy(ITK_Vol)


def cannyEdge(vol, var, maxErr, minT, maxT):
    ITK_Vol = imageFromNumpyToITK(vol)
    cannyEdge = sitk.CannyEdgeDetectionImageFilter()
    cannyEdge.SetLowerThreshold(minT)
    cannyEdge.SetUpperThreshold(maxT)
    cannyEdge.SetVariance([var, var, var])
    cannyEdge.SetMaximumError(maxErr)
    ITK_Vol = cannyEdge.Execute(ITK_Vol)
    return imageFromITKToNumpy(ITK_Vol)


def recursiveGauss(vol, sigma):
    ITK_Vol = imageFromNumpyToITK(vol)
    recurGaussX = sitk.RecursiveGaussianImageFilter()
    recurGaussY = sitk.RecursiveGaussianImageFilter()
    recurGaussX.SetSigma(sigma)
    recurGaussY.SetSigma(sigma)
    recurGaussY.SetDirection(1)
    ITK_Vol = recurGaussY.Execute(recurGaussX.Execute(ITK_Vol))

    return imageFromITKToNumpy(ITK_Vol)


def median(vol, radius, flag_2D=False):
    if flag_2D == 0:
        ITK_Vol = imageFromNumpyToITK(vol)
        med = sitk.MedianImageFilter()
        med.SetRadius(radius)
        ITK_Vol = med.Execute(ITK_Vol)
        Vol = imageFromITKToNumpy(ITK_Vol)
    else:

        Vol = np.zeros((vol.shape[0], vol.shape[1], vol.shape[2]))
        for z in range(0, vol.shape[0]):
            Vol[z, :, :] = ndimage.median_filter(vol[z, :, :], 3)

    return Vol


def sigmo(vol, alpha, beta):
    ITK_Vol = imageFromNumpyToITK(vol)
    return imageFromITKToNumpy(sitk.Sigmoid(ITK_Vol, alpha, beta))


def gradGauss(vol, sigma):
    ITK_Vol = imageFromNumpyToITK(vol)
    return imageFromITKToNumpy(sitk.GradientMagnitudeRecursiveGaussian(ITK_Vol, sigma))


def FFTCOR(ImF, ImM, sizeF, sizeM, overlap):
    ITK_M = imageFromNumpyToITK(ImM)
    ITK_F = imageFromNumpyToITK(ImF)

    mask_M = np.zeros((ImM.shape[0], ImM.shape[1], ImM.shape[2]))
    mask_F = np.zeros((ImF.shape[0], ImF.shape[1], ImF.shape[2]))

    mask_M[50:100, 50:100, 50:100] = 1.0
    mask_F[50:100, 50:100, 50:100] = 1.0

    ITK_mM = imageFromNumpyToITK(mask_M)
    ITK_mF = imageFromNumpyToITK(mask_F)

    filterCorr = sitk.MaskedFFTNormalizedCorrelationImageFilter()
    filterCorr.SetRequiredFractionOfOverlappingPixels(overlap)

    return imageFromITKToNumpy(filterCorr.Execute(ITK_F, ITK_M, ITK_mF, ITK_mM))


def morpho(Operation, vol, kernel):
    vol = np.copy(vol.astype(np.uint8))

    if Operation != 'Fill':
        vol_itk = imageFromNumpyToITK(vol)

        if Operation == 'Dilate':
            filterMorpho = sitk.BinaryDilateImageFilter()
        if Operation == 'Erode':
            filterMorpho = sitk.BinaryErodeImageFilter()
        if Operation == 'Open':
            filterMorpho = sitk.BinaryMorphologicalOpeningImageFilter()
        if Operation == 'Close':
            filterMorpho = sitk.BinaryMorphologicalClosingImageFilter()
        if Operation == 'Skeletton':
            filterMorpho = sitk.BinaryThinningImageFilter()

        if Operation != 'Skeletton':
            filterMorpho.SetKernelRadius(kernel)
        vol_out_itk = filterMorpho.Execute(vol_itk)
        return imageFromITKToNumpy(vol_out_itk)

    if Operation == 'Fill':

        shapeVol = vol.shape
        for z in range(0, shapeVol[0]):
            image = vol[z, :, :]
            imageITK = imageFromNumpyToITK(image)
            filterMorpho = sitk.BinaryFillholeImageFilter()
            imageITK = filterMorpho.Execute(imageITK)
            vol[z, :, :] = imageFromITKToNumpy(imageITK)

        return vol


def Distance(npArray):
    ITK_LS = imageFromNumpyToITK(npArray)
    distanceDetect = sitk.SignedDanielssonDistanceMapImageFilter()
    return imageFromITKToNumpy(distanceDetect.Execute(ITK_LS))


""" Utility """


def imageFromNumpyToITK(vol):
    return sitk.GetImageFromArray(vol)


def imageFromITKToNumpy(vol):
    return sitk.GetArrayFromImage(vol)


def returnZonesForEqualize(zones, data):
    if len(zones) == 1:
        x1 = zones[-1][0]
        x2 = zones[-1][3]
        y1 = zones[-1][1]
        y2 = zones[-1][4]
        z1 = zones[-1][2]
        z2 = zones[-1][5]
    elif len(zones) > 1:
        x1 = zones[-2][0]
        x2 = zones[-2][3]
        y1 = zones[-2][1]
        y2 = zones[-2][4]
        z1 = zones[-2][2]
        z2 = zones[-2][5]

    if z2 < z1:
        zs = z1
        z1 = z2
        z2 = zs
    if y2 < y1:
        ys = y1
        y1 = y2
        y2 = ys
    if x2 < x1:
        xs = x1
        x1 = x2
        x2 = xs

    if z1 < 0:
        z1 = 0
    elif z1 > data.shape[0]:
        msgBox = qt.QMessageBox()
        msgBox.setText('Out of Volume')
        msgBox.exec_()

    if z2 < 0:
        msgBox = qt.QMessageBox()
        msgBox.setText('Out of Volume')
        msgBox.exec_()
    elif z2 > data.shape[0]:
        z2 = data.shape[0] - 1

    if y1 < 0:
        y1 = 0
    elif y1 > data.shape[1]:
        msgBox = qt.QMessageBox()
        msgBox.setText('Out of Volume')
        msgBox.exec_()

    if y2 < 0:
        msgBox = qt.QMessageBox()
        msgBox.setText('Out of Volume')
        msgBox.exec_()
    elif y2 > data.shape[1]:
        y2 = data.shape[1] - 1

    if x1 < 0:
        x1 = 0
    elif x1 > data.shape[2]:
        msgBox = qt.QMessageBox()
        msgBox.setText('Out of Volume')
        msgBox.exec_()

    if x2 < 0:
        msgBox = qt.QMessageBox()
        msgBox.setText('Out of Volume')
        msgBox.exec_()
    elif x2 > data.shape[2]:
        x2 = data.shape[2] - 1

    if len(zones) == 1:
        zones_in = [[z1, y1, x1, z2, y2, x2]]

    elif len(zones) > 1:

        x12 = zones[-1][0]
        x22 = zones[-1][3]
        y12 = zones[-1][1]
        y22 = zones[-1][4]
        z12 = zones[-1][2]
        z22 = zones[-1][5]

        if z22 < z12:
            zs = z12
            z12 = z22
            z22 = zs
        if y22 < y12:
            ys = y12
            y12 = y22
            y22 = ys
        if x22 < x12:
            xs = x12
            x12 = x22
            x22 = xs

        if z12 < 0:
            z12 = 0
        elif z12 > data.shape[0]:
            msgBox = qt.QMessageBox()
            msgBox.setText('Out of Volume')
            msgBox.exec_()

        if z22 < 0:
            msgBox = qt.QMessageBox()
            msgBox.setText('Out of Volume')
            msgBox.exec_()
        elif z22 > data.shape[0]:
            z22 = data.shape[0] - 1

        if y12 < 0:
            y12 = 0
        elif y12 > data.shape[1]:
            msgBox = qt.QMessageBox()
            msgBox.setText('Out of Volume')
            msgBox.exec_()

        if y22 < 0:
            msgBox = qt.QMessageBox()
            msgBox.setText('Out of Volume')
            msgBox.exec_()
        elif y22 > data.shape[1]:
            y22 = data.shape[1] - 1

        if x12 < 0:
            x12 = 0
        elif x12 > data.shape[2]:
            msgBox = qt.QMessageBox()
            msgBox.setText('Out of Volume')
            msgBox.exec_()

        if x22 < 0:
            msgBox = qt.QMessageBox()
            msgBox.setText('Out of Volume')
            msgBox.exec_()
        elif x22 > data.shape[2]:
            x22 = data.shape[2] - 1

        zones_in = [[z1, y1, z1, z2, y2, x2], [z12, y12, z12, z22, y22, x22]]

    return zones_in


def saveImageAsTxt(I1, I2, I3, I4, I5, I6, filenames, binXls):
    sizeZ = I1.shape[0]
    sizeX = I2.shape[1]
    sizeY = I3.shape[2]

    with open(filenames[0], 'w'):
        pass
    with open(filenames[1], 'w'):
        pass
    with open(filenames[2], 'w'):
        pass
    with open(filenames[3], 'w'):
        pass
    with open(filenames[4], 'w'):
        pass
    with open(filenames[5], 'w'):
        pass
    with open(filenames[6], 'w'):
        pass
    with open(filenames[7], 'w'):
        pass
    with open(filenames[8], 'w'):
        pass

    data_to_save1 = np.empty((10000000, 1))
    data_to_save2 = np.empty((10000000, 1))
    data_to_save3 = np.empty((10000000, 1))
    data_to_save4 = np.empty((10000000, 1))
    data_to_save5 = np.empty((10000000, 1))
    data_to_save6 = np.empty((10000000, 1))
    data_to_save7 = np.empty((10000000, 1))
    data_to_save8 = np.empty((10000000, 1))
    data_to_save9 = np.empty((10000000, 1))

    index = 0

    for z in range(0, sizeZ, binXls):
        for x in range(0, sizeX, binXls):
            for y in range(0, sizeY, binXls):
                if I1[z, x, y] != -2000.0:
                    I1px = I1[z, x, y]
                    I2px = I2[z, x, y]
                    I3px = I3[z, x, y]
                    I4px = I4[z, x, y]
                    I5px = I5[z, x, y]
                    I6px = I6[z, x, y]

                    data_to_save1[index] = x
                    data_to_save2[index] = y
                    data_to_save3[index] = z
                    data_to_save4[index] = I1px
                    data_to_save5[index] = I2px
                    data_to_save6[index] = I3px
                    data_to_save7[index] = I4px
                    data_to_save8[index] = I5px
                    data_to_save9[index] = I6px

                    index += 1

    np.savetxt(filenames[0], data_to_save1[0:index], delimiter='\t')
    np.savetxt(filenames[1], data_to_save2[0:index], delimiter='\t')
    np.savetxt(filenames[2], data_to_save3[0:index], delimiter='\t')
    np.savetxt(filenames[3], data_to_save4[0:index], delimiter='\t')
    np.savetxt(filenames[4], data_to_save5[0:index], delimiter='\t')
    np.savetxt(filenames[5], data_to_save6[0:index], delimiter='\t')
    np.savetxt(filenames[6], data_to_save7[0:index], delimiter='\t')
    np.savetxt(filenames[7], data_to_save8[0:index], delimiter='\t')
    np.savetxt(filenames[8], data_to_save9[0:index], delimiter='\t')


def hessian_matrix(x):
    x_grad = np.gradient(x)
    hessian = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype)
    for k, grad_k in enumerate(x_grad):
        tmp_grad = np.gradient(grad_k)
        for l, grad_kl in enumerate(tmp_grad):
            hessian[k, l, :, :] = grad_kl
    return hessian


def computeEigenvalueHessianMatrix(H):
    i, j, x, y, z = H.shape
    H_reshape = H.reshape(i, j, x * y * z)
    H_reshape = np.moveaxis(H_reshape, -1, 0)
    length = H_reshape.shape[0]
    #Eig = np.apply_along_axis(LA.eigh,0,H_reshape.T)


    image_lambda_1 = np.zeros(H_reshape.shape[0])
    image_lambda_2 = np.zeros(H_reshape.shape[0])
    image_lambda_3 = np.zeros(H_reshape.shape[0])



    for i in range(H_reshape.shape[0]):
        v ,_ = LA.eigh(H_reshape[i])
        image_lambda_1[i] = v[0]
        image_lambda_2[i] = v[1]
        image_lambda_3[i] = v[2]

        print((i/length)*100.0)

    image_lambda_1 = image_lambda_1.reshape(x,y,z)
    image_lambda_2 = image_lambda_2.reshape(x, y, z)
    image_lambda_3 = image_lambda_3.reshape(x, y, z)

    return [image_lambda_1,image_lambda_2,image_lambda_3]


def power_method(A, iterations):
    b_k = np.random.rand(A.shape[1], A.shape[2])
    b_kn = b_k
    norm_vec = np.vectorize(np.linalg.norm)
    for _ in range(iterations):
        b_k1 = np.matmul(A, b_k)
        b_kn1 = np.dot(A[100],b_kn)
        print('--------------')
        print('Dot')
        print(b_k1[100])
        print(b_kn1)
        b_k = b_k1/(norm_vec(b_k1))
        b_kn = b_kn1/np.linalg.norm(b_kn1)
        print('-------------')
        print('Norm')
        print(b_k[100])
        print(b_kn)



    #b_ki = np.nan_to_num(b_ki)
    x = np.matmul(A, X0)
    norm = np.linalg.norm(x, axis=0)


    return x/norm


def frangi(nd_array, scale_range=(1, 10), scale_step=2, alpha=0.5, beta=0.5, frangi_c=500, black_vessels=True):

    if not nd_array.ndim == 3:
        raise(ValueError("Only 3 dimensions is currently supported"))

    # from https://github.com/scikit-image/scikit-image/blob/master/skimage/filters/_frangi.py#L74
    sigmas = np.arange(scale_range[0], scale_range[1], scale_step)
    if np.any(np.asarray(sigmas) < 0.0):
        raise ValueError("Sigma values less than zero are not valid")

    filtered_array = np.zeros(sigmas.shape + nd_array.shape)

    for i, sigma in enumerate(sigmas):
        print(i)
        eigenvalues = absolute_hessian_eigenvalues(nd_array, sigma=sigma, scale=True)
        filtered_array[i] = compute_vesselness(*eigenvalues, alpha=alpha, beta=beta, c=frangi_c,
                                               black_white=black_vessels)

    return np.max(filtered_array, axis=0)


def compute_measures(eigen1, eigen2, eigen3):
    """
    RA - plate-like structures
    RB - blob-like structures
    S - background
    """
    Ra = divide_nonzero(np.abs(eigen2), np.abs(eigen3))
    Rb = divide_nonzero(np.abs(eigen1), np.sqrt(np.abs(np.multiply(eigen2, eigen3))))
    S = np.sqrt(np.square(eigen1) + np.square(eigen2) + np.square(eigen3))
    return Ra, Rb, S


def compute_plate_like_factor(Ra, alpha):
    return 1 - np.exp(np.negative(np.square(Ra)) / (2 * np.square(alpha)))


def compute_blob_like_factor(Rb, beta):
    return np.exp(np.negative(np.square(Rb) / (2 * np.square(beta))))


def compute_background_factor(S, c):
    return 1 - np.exp(np.negative(np.square(S)) / (2 * np.square(c)))


def compute_vesselness(eigen1, eigen2, eigen3, alpha, beta, c, black_white):
    Ra, Rb, S = compute_measures(eigen1, eigen2, eigen3)
    plate = compute_plate_like_factor(Ra, alpha)
    blob = compute_blob_like_factor(Rb, beta)
    background = compute_background_factor(S, c)
    return filter_out_background(plate * blob * background, black_white, eigen2, eigen3)


def filter_out_background(voxel_data, black_white, eigen2, eigen3):
    """
    Set black_white to true if vessels are darker than the background and to false if
    vessels are brighter than the background.
    """
    if black_white:
        voxel_data[eigen2 < 0] = 0
        voxel_data[eigen3 < 0] = 0
    else:
        voxel_data[eigen2 > 0] = 0
        voxel_data[eigen3 > 0] = 0
    voxel_data[np.isnan(voxel_data)] = 0
    return voxel_data

def compute_hessian_matrix(nd_array, sigma=1, scale=True):
    """
    Computes the hessian matrix for an nd_array.
    This can be used to detect vesselness as well as other features.
    In 3D the first derivative will contain three directional gradients at each index:
    [ gx,  gy,  gz ]
    The Hessian matrix at each index will then be equal to the second derivative:
    [ gxx, gxy, gxz]
    [ gyx, gyy, gyz]
    [ gzx, gzy, gzz]
    The Hessian matrix is symmetrical, so gyx == gxy, gzx == gxz, and gyz == gzy.
    :param nd_array: n-dimensional array from which to compute the hessian matrix.
    :param sigma: gaussian smoothing to perform on the array.
    :param scale: if True, the hessian elements will be scaled by sigma squared.
    :return: hessian array of shape (..., ndim, ndim)
    """
    ndim = nd_array.ndim

    # smooth the nd_array
    smoothed = ndi.gaussian_filter(nd_array, sigma=sigma)

    # compute the first order gradients
    gradient_list = np.gradient(smoothed)

    # compute the hessian elements
    hessian_elements = [np.gradient(gradient_list[ax0], axis=ax1)
                        for ax0, ax1 in combinations_with_replacement(range(ndim), 2)]

    if sigma > 0 and scale:
        # scale the elements of the hessian matrix
        hessian_elements = [(sigma ** 2) * element for element in hessian_elements]

    # create hessian matrix from hessian elements
    hessian_full = [[None] * ndim] * ndim

    for index, (ax0, ax1) in enumerate(combinations_with_replacement(range(ndim), 2)):
        element = hessian_elements[index]
        hessian_full[ax0][ax1] = element
        if ax0 != ax1:
            hessian_full[ax1][ax0] = element

    hessian_rows = list()
    for row in hessian_full:
        hessian_rows.append(np.stack(row, axis=-1))

    hessian = np.stack(hessian_rows, axis=-2)
    return hessian


def absolute_hessian_eigenvalues(nd_array, sigma=1, scale=True):
    """
    Eigenvalues of the hessian matrix calculated from the input array sorted by absolute value.
    :param nd_array: input array from which to calculate hessian eigenvalues.
    :param sigma: gaussian smoothing parameter.
    :param scale: if True hessian values will be scaled according to sigma squared.
    :return: list of eigenvalues [eigenvalue1, eigenvalue2, ...]
    """
    return absolute_eigenvaluesh(compute_hessian_matrix(nd_array, sigma=sigma, scale=scale))
if __name__ == "__main__":

    Image = np.random.rand(150,150,150)
    filtered = frangi(Image)
    print(filtered.shape)
