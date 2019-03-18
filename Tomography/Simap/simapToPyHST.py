# -----------------------------------------------------------------------
# Copyright: 2010-2018, imec Vision Lab, University of Antwerp
#            2013-2018, CWI, Amsterdam
#
# Contact: astra@astra-toolbox.com
# Website: http://www.astra-toolbox.com/
#
# This file is part of the ASTRA Toolbox.
#
#
# The ASTRA Toolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ASTRA Toolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.
#
# -----------------------------------------------------------------------

import numpy as np
from numpy import linalg as LA
from skimage import io
import pylab
import os


def loadgeom(filename):
    g = np.loadtxt(filename, dtype={
        'names': ('name', 'Sx', 'Sy', 'Sz', 'P1x', 'P1y', 'P1z', 'P2x', 'P2y', 'P2z', 'P3x', 'P3y', 'P3z'),
        'formats': ('U80', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f')}, delimiter=';', skiprows=2)
    geo = {}
    for p in g:
        l = list(p)
        geo[p[0]] = l[1:]
    return geo


class scan:

    def __init__(self, param):
        self.direc = param['direc']
        self.prefix = param['prefix']
        self.direcout = param['direcout']
        self.prefixout = param['prefixout']
        self.detector = param['detector']
        self.projs = param['projs']
        self.volume = param['volume']
        self.resolution = param['resolution']
        self.correction = param['correction']
        self.algo = param['algo']
        self.vec = param['vec']
        self.show = param['show']

    def convertgeom(self, geom):
        nbproj = len(self.projs['views'])
        if (self.vec):
            if len(self.projs['frames']) > 0:
                l = geom["%s_%04d_%04d.tiff" % (self.prefix, self.projs['views'][0], 0)]
            else:
                l = geom["%s%05d.tif" % (self.prefix, self.projs['views'][0])]

            SP1 = np.array(l[3:6]) - np.array(l[0:3])

            P1P3 = np.array(l[6:9]) - np.array(l[3:6])
            nP1P3 = LA.norm(P1P3)
            u = P1P3 / nP1P3
            SOURCE_X = np.inner(SP1, u) * self.detector['col_count'] / (2. * nP1P3) + (
                        self.detector['col_count'] - 1) / 2.

            P1P2 = np.array(l[9:12]) - np.array(l[3:6])
            nP1P2 = LA.norm(P1P2)
            v = P1P2 / nP1P2
            SOURCE_Y = np.inner(SP1, v) * self.detector['row_count'] / (2. * nP1P2) + (
                        self.detector['row_count'] - 1) / 2.

            SP0 = SP1 - np.inner(SP1, u) * u - np.inner(SP1, v) * v

            ez = SP0 / LA.norm(SP0)
            ey = np.array((0, 1, 0))
            ex = np.cross(ey, ez)

            SO = -np.array(l[0:3])
            OFFSET_X = np.inner(SO, ex)
            OFFSET_Y0 = np.inner(SO, ey)

            sdd = LA.norm(SP0)
            SOURCE_DISTANCE = np.inner(SO, ez)
            DETECTOR_DISTANCE = sdd - SOURCE_DISTANCE

            angles = np.zeros((len(self.projs['views']),))
            axis_corrections = np.zeros((len(self.projs['views']), 2))
            for i, p in zip(np.arange(0, len(self.projs['views'])), self.projs['views']):
                if len(self.projs['frames']) > 0:
                    l = geom["%s_%04d_%04d.tiff" % (self.prefix, p, 0)]
                else:
                    l = geom["%s%05d.tif" % (self.prefix, p)]

                SP1 = np.array(l[3:6]) - np.array(l[0:3])

                P1P3 = np.array(l[6:9]) - np.array(l[3:6])
                nP1P3 = LA.norm(P1P3)
                u = P1P3 / nP1P3

                P1P2 = np.array(l[9:12]) - np.array(l[3:6])
                nP1P2 = LA.norm(P1P2)
                v = P1P2 / nP1P2

                SP0 = SP1 - np.inner(SP1, u) * u - np.inner(SP1, v) * v

                ez = SP0 / LA.norm(SP0)
                ey = np.array((0, 1, 0))
                ex = np.cross(ey, ez)

                SO = -np.array(l[0:3])
                OFFSET_X = np.inner(SO, ex)
                OFFSET_Y = np.inner(SO, ey) - OFFSET_Y0
                axis_corrections[i, :] = np.array((OFFSET_X, OFFSET_Y))
                sdd = LA.norm(SP0)
                sod = np.inner(SO, ex)
                odd = sdd - sod

                angles[i] = np.arccos(ez[0]) / np.pi * 180. * np.sign(ex[0])
                # dy=(odd*astra_geom[i,2]+sod*astra_geom[i,5])/(sod+odd)+self.detector['offsetImY']*self.detector['pixsize']/(1.+odd/sod)

            pyhst2_geom = {'SOURCE_DISTANCE': SOURCE_DISTANCE / 1000, \
                           'DETECTOR_DISTANCE': DETECTOR_DISTANCE / 1000, \
                           'SOURCE_X': SOURCE_X, \
                           'SOURCE_Y': SOURCE_Y, \
                           'AXIS_CORRECTIONS': axis_corrections, \
                           'ANGLES': angles}

            print(pyhst2_geom)

            pyhst2_geom_corr = self.correctgeom(pyhst2_geom)

            if (self.show['Geom']):
                print(pyhst2_geom)
            return pyhst2_geom
        else:
            angles = np.zeros((len(self.projs['views']),))
            l = geom["%s_%04d_%04d.tiff" % (self.prefix, self.projs['views'][0], 0)]
            sod = np.sqrt(l[0] ** 2 + l[2] ** 2)
            odd = np.sqrt(l[3] ** 2 + l[5] ** 2)
            for i, p in zip(np.arange(0, len(self.projs['views'])), self.projs['views']):
                l = geom["%s_%04d_%04d.tiff" % (self.prefix, p, 0)]
                angles[i] = np.arccos(l[2] / sod) * np.sign(l[0])
            if (self.show['Geom']):
                print(sod, odd, angles / 3.14 * 180)
            return pyhst_geom

    def correctgeom(self, geom):
        # geom[:,0]-=self.correction['offsetX']*self.detector['pixsize']
        # geom[:,3]-=self.correction['offsetX']*self.detector['pixsize']
        return geom

    def getgeom(self):
        return self.convertgeom(loadgeom("%s/%s.csv" % (self.direc, self.prefix)))

    def getproj(self):
        direc = self.direc
        prefix = self.prefix
        projs = self.projs['views']
        frames = self.projs['frames']
        detector = self.detector

        for i, p in zip(np.arange(0, len(projs)), projs):
            if (len(frames) > 1):
                fp = np.zeros((len(frames), detector['roiY'], detector['roiX']))
                for j, f in zip(np.arange(0, len(frames)), frames):
                    fpf = io.imread("%s/%s_%04d_%04d.tiff" % (direc, prefix, p, f))
                    fp[j, :, :] = fpf[int(detector['row_count'] / 2 - detector['roiY'] / 2 - detector['offsetImY']): \
                                      int(detector['row_count'] / 2 + detector['roiY'] / 2 - detector['offsetImY']), \
                                  int(detector['col_count'] / 2 - detector['roiX'] / 2 - detector['offsetImX']): \
                                  int(detector['col_count'] / 2 + detector['roiX'] / 2 - detector['offsetImX'])]
                if (self.projs['fusion'] == 'med'):
                    proj = np.median(fp, axis=0)
                else:
                    proj = np.mean(fp, axis=0)
            else:
                if (len(frames) > 0):
                    fpf = io.imread("%s/%s_%04d_%04d.tiff" % (direc, prefix, p, 0))
                else:
                    fpf = io.imread("%s/%s%05d.tif" % (direc, prefix, p))

                proj = fpf[int(detector['row_count'] / 2 - detector['roiY'] / 2 - detector['offsetImY']): \
                           int(detector['row_count'] / 2 + detector['roiY'] / 2 - detector['offsetImY']), \
                       int(detector['col_count'] / 2 - detector['roiX'] / 2 - detector['offsetImX']): \
                       int(detector['col_count'] / 2 + detector['roiX'] / 2 - detector['offsetImX'])]

            # np.float32(proj)

            newedf = EdfFile("%s/%s_%04d.edf" % (self.direcout, self.prefixout, i))
            newedf.WriteImage({}, proj, Append=0, DataType='UnsignedShort', ByteOrder="LowByteFirst")
            del newedf
            if i in self.show['Proj']:
                pylab.figure()
                pylab.imshow(proj)
        return "%s/%s_" % (self.direcout, self.prefixout)

    # def correctproj(self, proj):
    #    proj_corr=proj
    #    imax = np.max(np.max(proj))
    #    imax = 250

    #   if (self.correction['log']):
    #        proj_corr = -np.log((proj+0.01)/imax)

    #    for i in self.show['Projcor']:
    #        pylab.figure()
    #        pylab.imshow(proj_corr[:,i,:])

    #    return proj_corr

    def createANGLES(self, proj_geom, anglesname=""):
        if (anglesname == ""):
            anglesname = "%s/%s_angles.txt" % (self.direcout, self.prefixout)
        np.savetxt(anglesname, proj_geom['ANGLES'])

    def createAXISCORRECTIONS(self, proj_geom, axiscorrectionsname=""):
        if (axiscorrectionsname == ""):
            axiscorrectionsname = "%s/%s_axiscorrections.txt" % (self.direcout, self.prefixout)
        np.savetxt(axiscorrectionsname, proj_geom['AXIS_CORRECTIONS'], )

    def createPAR(self, proj_geom, projections_prefix, vol_prefix, parname, anglesname="", axiscorrectionsname=""):
        if (anglesname == ""):
            anglesname = "%s/%s_angles.txt" % (self.direcout, self.prefixout)
        if (axiscorrectionsname == ""):
            axiscorrectionsname = "%s/%s_axiscorrections.txt" % (self.direcout, self.prefixout)

        f = open(parname, 'w')

        f.write("# HST_SLAVE PARAMETER FILE\n\n")

        f.write("# Parameters defining the projection file series\n")
        f.write("FILE_PREFIX = %s\n" % (projections_prefix))
        f.write("NUM_FIRST_IMAGE = 0 # No. of first projection file\n")
        f.write("NUM_LAST_IMAGE = %d  # No. of last projection file\n" % (len(self.projs['views']) - 1))
        f.write("NUMBER_LENGTH_VARIES = NO\n")
        f.write("LENGTH_OF_NUMERICAL_PART = 4 # No. of characters\n")
        f.write("FILE_POSTFIX = .edf\n")
        f.write("FILE_INTERVAL = 1 # Interval between input files\n\n")

        f.write("# Parameters defining the projection file format\n")
        f.write("NUM_IMAGE_1 = %d  # Number of pixels horizontally\n" % (self.detector['roiX']))
        f.write("NUM_IMAGE_2 = %d  # Number of pixels vertically\n" % (self.detector['roiY']))
        f.write("IMAGE_PIXEL_SIZE_1 = %f # Pixel size horizontally (microns)\n" % (self.detector['pixsize'] * 1000.))
        f.write("IMAGE_PIXEL_SIZE_2 = %f # Pixel size vertically\n\n" % (self.detector['pixsize'] * 1000.))

        f.write("# Parameters defining background treatment\n")
        f.write("SUBTRACT_BACKGROUND = NO  # Subtract background from data\n")
        f.write("BACKGROUND_FILE = mydarks.edf\n\n")

        f.write("# Parameters defining flat-field treatment\n")
        f.write("CORRECT_FLATFIELD = NO  # Divide by flat-field image\n")
        f.write("FLATFIELD_CHANGING = NO  # Series of flat-field files\n")
        f.write("FLATFIELD_FILE = N.A.\n")
        f.write("FF_PREFIX = myflat\n")
        f.write("FF_NUM_FIRST_IMAGE = 0 # No. of first flat-field file\n")
        f.write("FF_NUM_LAST_IMAGE = 1500 # No. of last flat-field file\n")
        f.write("FF_NUMBER_LENGTH_VARIES = NO\n")
        f.write("FF_LENGTH_OF_NUMERICAL_PART = 4 # No. of characters\n")
        f.write("FF_POSTFIX = .edf\n")
        f.write("FF_FILE_INTERVAL = 1500 # Interval between flat-field files\n\n")

        f.write("ZEROCLIPVALUE = 0.01 # Minimum value of radiographs after flat / before log\n")
        f.write("#ONECLIPVALUE\n\n")

        f.write("# Parameters defining experiment\n")
        f.write("ANGLES_FILE = %s\n" % (anglesname))
        f.write("ANGLE_BETWEEN_PROJECTIONS = 10 # Increment angle in degrees\n")
        f.write("ROTATION_VERTICAL = YES\n")
        f.write("ROTATION_AXIS_POSITION = %f       # Position in pixels\n\n" % (
                    self.detector['roiX'] / 2. + self.detector['offsetImX'] + self.correction['offsetX']))

        f.write("# Parameters defining reconstruction\n")
        f.write("OUTPUT_SINOGRAMS = NO # Output sinograms to files or not\n")
        f.write("OUTPUT_RECONSTRUCTION = YES # Reconstruct and save or not\n")
        f.write("START_VOXEL_1 =      %d # X-start of reconstruction volume\n" % (self.volume['offsetX']))
        f.write("START_VOXEL_2 =      %d # Y-start of reconstruction volume\n" % (self.volume['offsetY']))
        f.write("START_VOXEL_3 =     %d # Z-start of reconstruction volume\n" % (self.volume['offsetZ']))
        f.write("END_VOXEL_1 =    %d  # X-end of reconstruction volume\n" % (
                    self.volume['offsetX'] + self.volume['dimx'] - 1))
        f.write("END_VOXEL_2 =    %d  # Y-end of reconstruction volume\n" % (
                    self.volume['offsetY'] + self.volume['dimy'] - 1))
        f.write("END_VOXEL_3 =    %d  # Z-end of reconstruction volume\n" % (
                    self.volume['offsetZ'] + self.volume['dimz'] - 1))

        f.write("OVERSAMPLING_FACTOR = 4 # 0 = Linear, 1 = Nearest pixel\n")
        f.write("ANGLE_OFFSET = 0.000000 # Reconstruction rotation offset angle in degrees\n")
        f.write("CACHE_KILOBYTES = 1024 # Size of processor cache (L2) per processor (KBytes)\n")
        f.write("SINOGRAM_MEGABYTES = 2000 # Maximum size of sinogram storage (megabytes)\n\n")

        f.write("# Parameters extra features PyHST\n")
        f.write("DO_CCD_FILTER = 1 # CCD filter (spikes)\n")
        f.write('CCD_FILTER = "CCD_Filter"\n')
        f.write('CCD_FILTER_PARA = {"threshold": 0.040000 }\n')
        f.write("DO_SINO_FILTER = 1 # Sinogram filter (rings)\n")
        f.write('SINO_FILTER = "SINO_Filter"\n')
        f.write("ar = numpy.ones(2048,'f')\n")
        f.write("ar[0]=0.0\n")
        f.write("ar[2:18]=0.0\n")
        f.write('SINO_FILTER_PARA = {"FILTER": ar }\n\n')

        f.write("#DISTORSION_FILE = None\n")
        f.write(
            "#The horizontal distortion will be read from the file mydist_H.edf and the vertical one from mydist_V.edf\n\n")

        f.write("DO_AXIS_CORRECTION = YES # Axis correction\n")
        f.write("AXIS_CORRECTION_FILE = %s #in pixel positive in ccd axis directions\n" % (axiscorrectionsname))
        f.write('''OPTIONS= { 'padding':'E' , 'axis_to_the_center':'Y'} # Padding and position axis\n\n''')

        f.write("#BH_LUT_FILE = None # 2 columns : integral of mu and multiplicative factor after flat field\n\n")
        f.write("")
        f.write("# Parameters for Paganin reconstruction\n")
        f.write("DO_PAGANIN = 0\n")
        f.write("PAGANIN_Lmicron = 247.790366 \n")
        f.write("PAGANIN_MARGE = 20 \n")
        f.write("DO_OUTPUT_PAGANIN = 0\n")
        f.write("OUTPUT_PAGANIN_FILE = projes/paga_cufft_\n")
        f.write("PAGANIN_TRY_CUFFT = 1\n")
        f.write("AGANIN_TRY_FFTW = 1\n\n")

        f.write("CONICITY = 1\n")
        f.write("CONICITY_FAN = 0\n")
        f.write("DETECTOR_DISTANCE = %f\n" % (proj_geom['DETECTOR_DISTANCE']))
        f.write("SOURCE_DISTANCE = %f\n" % (proj_geom['SOURCE_DISTANCE']))
        f.write("SOURCE_X = %f\n" % (proj_geom['SOURCE_X']))
        f.write("SOURCE_Y = %f\n\n" % (proj_geom['SOURCE_Y']))

        f.write("DXPERPROJ = 0.0 #helical\n")
        f.write("DZPERPROJ = 0.0 #helical\n\n")

        f.write("TAKE_LOGARITHM = NO  # Take log of projection values\n\n")

        f.write("# Parameters defining output file / format\n")
        f.write("OUTPUT_FILE =  %s.vol\n\n" % (vol_prefix))

        f.write("# Reconstruction program options\n")
        f.write("DISPLAY_GRAPHICS = NO # No images\n\n")

        f.write("PUS=0.0\n")
        f.write("PUC=0.0\n")
        f.write("UNSHARP_LoG=0\n\n")

        f.write("TRYGPU=1\n\n")

        f.write("#FBFILTER=0        # 0 = rampe \n\n")

        f.close()

    def runPyHST2(self, parname):
        os.system("%s %s" % (self.algo['pyhst'], parname))

    def recons(self):

        # Initialize pygpu
        # ctx = pygpu.init('cuda')
        # pygpu.set_default_context(ctx)

        # vol_geom = astra.create_vol_geom(self.volume['dimy'], self.volume['dimx'], self.volume['dimz'])

        # vol_gpuarr = pygpu.gpuarray.zeros(astra.geom_size(vol_geom), dtype='float32')

        # z, y, x = vol_gpuarr.shape
        # vol_link = astra.data3d.GPULink(vol_gpuarr.gpudata, x, y, z,vol_gpuarr.strides[-2])

        # Create a data object for the reconstruction
        # rec_id = astra.data3d.create('-vol', vol_geom)
        # rec_id = astra.data3d.link('-vol', vol_geom, vol_link)

        print("Get geometry")
        proj_geom = self.getgeom()

        print("Get projections")
        projections_name = self.getproj()

        print("Create angles file")
        self.createANGLES(proj_geom)

        print("Create axis correction file")
        self.createAXISCORRECTIONS(proj_geom)

        print("Create par file")
        par_name = "%s.par" % (self.prefixout)
        vol_name = "%s/%s" % (self.direcout, self.prefixout)
        self.createPAR(proj_geom, projections_name, vol_name, par_name)

        # residual_error = np.zeros(self.algo['nIters'])
        # for i in range(self.algo['nIters']):
        #    # Run a single iteration
        #    astra.algorithm.run(alg_id, 1)
        #    residual_error[i] = astra.algorithm.get_res_norm(alg_id)
        #    #rec = astra.data3d.get(rec_id)
        #    if i in self.show['Iter']:
        #        pylab.figure()
        #        pylab.imshow(vol_gpuarr[self.show['Slice'][0],100:-100,100:-100])
        # Get the result
        # rec = astra.data3d.get(rec_id)
        # for i in self.show['Slice']:
        #    pylab.figure()
        #    pylab.imshow(rec[i,100:-100,100:-100])

        # if(self.show['Residu']):
        #    pylab.figure()
        #    pylab.plot(residual_error)

        # Clean up. Note that GPU memory is tied up in the algorithm object,
        # and main RAM in the data objects.
        # astra.algorithm.delete(alg_id)
        # astra.data3d.delete(rec_id)
        # astra.data3d.delete(proj_id)

        print("Run PyHST2")
        self.runPyHST2(par_name)

        # self.showSlices()


def defaultscan():
    detector = {'name': 'Varian', 'row_count': 728, 'col_count': 920, 'pixsize': 0.254, 'roiX': 920, 'roiY': 728,
                'offsetImX': 0, 'offsetImY': 0}
    projs = {'views': np.arange(0, 35), 'frames': (0, 1, 2), 'fusion': 'med'}
    direc = 'testdata'
    direcout = 'result'
    prefix = 'tot2_param0001'
    prefixout = prefix
    volume = {'dimx': 800, 'dimy': 800, 'dimz': 8, 'offsetX': 0, 'offsetY': 0, 'offsetZ': 524}
    resolution = 0.02
    algo = {'pyhst': 'PyHST2_2018a', 'algo': 'SIRT3D_CUDA', 'nIters': 10, 'MinConstraint': 0., 'MaxConstraint': 0.03}
    correction = {'log': True, 'offsetX': 1.5}
    vec = True
    show = {'Proj': (0,), 'Projcor': (0,), 'Slice': (0,), 'Residu': True, 'Geom': False,
            'Iter': (0, 1, 2, 10, 20, 30, 40, 49)}

    param = {'direc': direc, 'prefix': prefix, 'direcout': direcout, 'prefixout': prefixout, 'detector': detector,
             'projs': projs, 'volume': volume, 'resolution': resolution, 'correction': correction, 'algo': algo,
             'vec': vec, 'show': show}
    return scan(param)


def marcscan():
    detector = {'name': 'Varian', 'row_count': 1456, 'col_count': 1840, 'pixsize': 0.127, 'roiX': 1840, 'roiY': 1456,
                'offsetImX': 0, 'offsetImY': 0}
    projs = {'views': np.linspace(0, 2399, 10), 'frames': (), 'fusion': 'med'}
    direc = '../../2018_11_08_Marc_good/Proj'
    direcout = '../../2018_11_08_Marc_good/result'
    prefix = '2018_11_08_Marc_good'
    prefixout = prefix
    volume = {'dimx': 1800, 'dimy': 1800, 'dimz': 8, 'offsetX': 0, 'offsetY': 0, 'offsetZ': 1000}
    resolution = 0.02
    algo = {'pyhst': 'PyHST2_2018a', 'algo': 'SIRT3D_CUDA', 'nIters': 10, 'MinConstraint': 0., 'MaxConstraint': 0.03}
    correction = {'log': True, 'offsetX': 0}
    vec = True
    show = {'Proj': (0,), 'Projcor': (0,), 'Slice': (0,), 'Residu': True, 'Geom': False,
            'Iter': (0, 1, 2, 10, 20, 30, 40, 49)}

    param = {'direc': direc, 'prefix': prefix, 'direcout': direcout, 'prefixout': prefixout, 'detector': detector,
             'projs': projs, 'volume': volume, 'resolution': resolution, 'correction': correction, 'algo': algo,
             'vec': vec, 'show': show}
    return scan(param)


def run():
    # detector = {'name':'Varian', 'row_count':1840, 'col_count':1560, 'pixsize':0.127}
    # detector = {'name':'Varian', 'row_count':728, 'col_count':920, 'pixsize':0.254}
    # s=defaultscan()
    s = marcscan()
    s.recons()

    pylab.show()


if __name__ == '__main__':
    loadgeom('/')
    # def main():
    import sys
#    run()

# main()

# sys.exit()

