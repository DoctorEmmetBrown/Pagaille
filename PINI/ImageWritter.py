# -*- coding: utf-8 -*-


from EdfFile import EdfFile
from TiffI0 import TiffIO
import numpy as np
import scipy.io as sio
import scipy
import SimpleITK as sitk

import datetime, time
import nibabel as nib


class ImageWritter(object) : 

    def __init__(self,filename,data) : 
        
        self.fileName = filename

        self.data=data
        self.slice=None

    def writeData(self) : 
            
        if self.fileName.endswith('.edf')  : 
            fileEdf = EdfFile(self.fileName, access='wb')
            fileEdf.WriteImage({}, self.data)

        if self.fileName.endswith('.mat')  : 
            NameFile = self.fileName.split('/')[-1]
            NameFile = NameFile.split('.')[0]
            sio.savemat(self.fileName,{NameFile:self.data})

        if self.fileName.endswith('.tiff')  : 
            tifImage = TiffIO.TiffIO(self.fileName, 'wb+')
            tifImage.writeImage(self.data)
   
        if self.fileName.endswith('.png')  : 
            scipy.misc.imsave(self.fileName, self.data)
            
        if self.fileName.endswith('.dcm')  : 
            write_dicom(self.data, self.fileName)
            
        if self.fileName.endswith('.npy')  : 
            np.save( self.fileName, self.data)



def write_nifti( pixel_array, filename):
    new_image = nib.Nifti1Image(pixel_array, affine=np.eye(4))
    nib.save(new_image, filename)
    

def write_dicom(pixel_array,filename):
    writer = sitk.ImageFileWriter()
    # Use the study/seriers/frame of reference information given in the meta-data
    # dictionary and not the automatically generated information from the file IO
    image_slice = sitk.GetImageFromArray(pixel_array)

    writer.KeepOriginalImageUIDOn()
    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")

    # Copy the meta-data except the rescale-intercept, rescale-slope

    # Set relevant keys indicating the change, modify or remove private tags as needed
    
    
    image_slice.SetMetaData("0008|0031", modification_time)
    image_slice.SetMetaData("0008|0021", modification_date)
    image_slice.SetMetaData("0008|0008", "DERIVED\SECONDARY")
    # Each of the UID components is a number (cannot start with zero) and separated by a '.'
    # We create a unique series ID using the date and time.
    image_slice.SetMetaData("0020|000e", "1.2.826.0.1.3680043.2.1125."+modification_date+".1"+modification_time)
    # Write to the output directory and add the extension dcm if not there, to force writing is in DICOM format.
    writer.SetFileName(filename)
    writer.Execute(image_slice)
    return




















































