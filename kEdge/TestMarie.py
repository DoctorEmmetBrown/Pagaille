import EdfFile
import numpy as np
from scipy.signal import convolve2d
import glob
import os
#from skimage.restoration import inpaint
import scipy
from scipy import ndimage as nd

def createNans(shape, dtype=float):
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a


def inpaint_nans(im):
    ipn_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]]) # kernel for inpaint_nans
    nans = np.isnan(im)
    while np.sum(nans)>0:
        print 'ici'
        im[nans] = 0
        vNeighbors = convolve2d((nans==False),ipn_kernel,mode='same',boundary='symm')
        im2 = convolve2d(im,ipn_kernel,mode='same',boundary='symm')
        im2[vNeighbors>0] = im2[vNeighbors>0]/vNeighbors[vNeighbors>0]
        im2[vNeighbors==0] = np.nan
        im2[(nans==False)] = im[(nans==False)]
        im = im2
        nans = np.isnan(im)
    return im



def fill_inpaint(array,invalid=None,max_iter=5,tol=0.5,kernel_size=1,method='localmean'):
    """Replace NaN elements in an array using an iterative image inpainting algorithm.
    The algorithm is the following:
    1) For each element in the input array, replace it by a weighted average
    of the neighbouring elements which are not NaN themselves. The weights depends
    of the method type. If ``method=localmean`` weight are equal to 1/( (2*kernel_size+1)**2 -1 )
    2) Several iterations are needed if there are adjacent NaN elements.
    If this is the case, information is "spread" from the edges of the missing
    regions iteratively, until the variation is below a certain threshold.

    Parameters
    ----------
    array : 2d np.ndarray
    invalid : mask for array (masked elements True), optional
    an array containing NaN elements that have to be replaced
    max_iter : int
    the number of iterations
    kernel_size : int
    the size of the kernel, default is 1
    method : str
    the method used to replace invalid values. Valid options are
    `localmean`, 'idw'.
    Returns
    -------
    filled : 2d np.ndarray
    a copy of the input array, where NaN elements have been replaced.
    Credits
    -------
    http://stackoverflow.com/a/17125125/512111
    """
    array = np.asarray(array,np.float64)

    if invalid is not None:
        array[np.where(invalid)] = np.nan

    # convert masked array to array filled with nans
    array = np.ma.filled(array,np.nan)

    filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64 )

    # indices where array is NaN
    inans, jnans = np.nonzero(np.isnan(array))

    # number of NaN elements
    n_nans = len(inans)

    # arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=np.float64)
    replaced_old = np.zeros( n_nans, dtype=np.float64)

    # depending on kernel type, fill kernel array
    if method == 'localmean':

        print 'kernel_size', kernel_size
        for i in xrange(2*kernel_size+1):
            for j in xrange(2*kernel_size+1):
                kernel[i,j] = 1
        print(kernel, 'kernel')

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
            [0.5,0.75,0.75,0.75,0.5],
            [0.5,0.75,1,0.75,0.5],
            [0.5,0.75,0.75,0.5,1],
            [0, 0.5, 0.5 ,0.5 ,0]])
        print(kernel, 'kernel')
    else:
        raise NotImplementedError('Method not valid. Should be one of [\'localmean\',\'idw\'].')

    # fill new array with input elements
    for i in xrange(array.shape[0]):
        for j in xrange(array.shape[1]):
            filled[i,j] = array[i,j]

    # make several passes
    # until we reach convergence
    for it in xrange(max_iter):
        print 'iteration', it
        # for each NaN element
        for k in xrange(n_nans):
            i = inans[k]
            j = jnans[k]

            # initialize to zero
            filled[i,j] = 0.0
            n = 0

            # loop over the kernel
            for I in xrange(2*kernel_size+1):
                for J in xrange(2*kernel_size+1):

                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:

                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :

                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:

                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1*kernel[I,J]

            # divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = np.nan

        # check if mean square difference between values of replaced
        #elements is below a certain tolerance
        print('tolerance', np.mean( (replaced_new-replaced_old)**2 ))
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in xrange(n_nans):
                replaced_old[l] = replaced_new[l]

    return filled

if __name__ == "__main__":
    folder='/Volumes/VISITOR/md950/id17/pinosh/s7_/'
    outputFolder='/Volumes/VISITOR/md950/id17/pinosh/s8_Manu/'
    listOfFiles=glob.glob(folder+'*.edf')
    for fichier in listOfFiles:
        edf=EdfFile.EdfFile(fichier,access='r')
        image=edf.GetData(0)
        typeImage=image.dtype
        nbLines,nbColumn=image.shape
        newImage=createNans((nbLines,1296),np.float)
        newImage[:,0:255]=image[:,0:255]
        newImage[:,260:515]=image[:,256:511]
        newImage[:,520:775]=image[:,512:767]
        newImage[:,780:1035]=image[:,768:1023]
        newImage[:,1040:1295]=image[:,1024:1279]

        mask=np.isnan(newImage)

        newImage=fill_inpaint(newImage)

        #newImage=inpaint_nans(newImage)
        #newImage = inpaint.inpaint_biharmonic(newImage, mask, multichannel=True)

        newImage=np.asarray(newImage,np.uint16)
        outputFileName=outputFolder+'/modified_'+os.path.basename(fichier)
        filetoWrite = EdfFile.EdfFile(outputFileName, access='wb+')
        filetoWrite.WriteImage({}, newImage)






