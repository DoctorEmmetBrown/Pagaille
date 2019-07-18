'''
Created 09/01/2018 by Luca Fardin.
Based on the code by L. Broche
'''


import RegistrationLuca as Registration

import InputOutput_fabio as InputOutput
import MaskCreator
import numpy as np
import sys
import SimpleITK as sitk

path_fixed =  '/data/id17/broncho/fardin/ForInes/AboveSmallImages/'
path_moving = '/data/id17/broncho/fardin/ForInes/BelowSmallImages/'
#path_maskf = '/data/id17/broncho/fardin/Mask1/'
#path_maskm = '/data/id17/broncho/fardin/Mask2/'
output_Transform = '/data/id17/broncho/fardin/ForInes/TestTransform/'


#Load Images
IO=InputOutput.InputOutput()
image_fixed  = IO.ReadTifVolume(path_fixed)
image_moving = IO.ReadTifVolume(path_moving)
#InputMaskF   = IO.ReadEdfVolume(path_maskf)
#InputMaskM   = IO.ReadEdfVolume(path_maskm)

image_moving_m=np.copy(image_moving)
image_fixed_m=np.copy(image_fixed)

#image_fixed_m[image_fixed<40000]=0
#image_moving_m[image_moving<40000]=0

#image_fixed=image_fixed_f[400:1400+400,300:1400+300]
#image_moving=image_moving_f[400:1400+400,300:1400+300]

#Apply rigid transform to align the samples
reg_method0 = Registration.Registering()
reg_method0.load_image(image_fixed_m.astype(np.double),image_moving_m.astype(np.double))
#reg_method.initMask(InputMaskF,InputMaskM)
reg_method0.init_Translation()
reg_method0.init_metric("Mattes Mutual Information",[50],True,0.5) #100
reg_method0.init_Interpolator("Linear Interpolation")
reg_method0.initOptimizer("LBFGSB",[1e-7,500,20,2000,1e+1])
reg_method0.initScaling([1],[0])
sys.stdout.flush()
reg_method0.Execute()
#reg_method.SaveTransform(output_Transform)
image_moving=reg_method0.Apply_Step(image_moving)
transform0=reg_method0.GetTransform()
del reg_method0


#OOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooo
#This part is to concatenate a second transform.
#There should be a smarter way to do that. In case we can try

reg_method = Registration.Registering()
reg_method.load_image(image_fixed,image_moving)
#reg_method.initMask(InputMaskF,InputMaskM)
sys.stdout.flush()
reg_method.init_metric("Means Squares",[],True,0.5) #100
reg_method.init_Interpolator("Linear Interpolation")
#reg_method.SetInitialTransform(transform0)
reg_method.init_Euler()
reg_method.initOptimizer("LBFGSB",[1e-8,500,20,2000,1e+0])
reg_method.initScaling([1],[0])
sys.stdout.flush()
#reg_method.Execute()
#reg_method.SaveTransform(output_Transform)
#image_moving=reg_method.Apply_Step()
#transform1=reg_method.GetTransform()
del reg_method



composite_transform=sitk.Transform(transform0)
#composite_transform.AddTransform(transform1)
#composite_transform.AddTransform(transform0)

#OOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooo

#Here I save the transform and the images.
sitk.WriteTransform(composite_transform,output_Transform+"Transform.tfm")
IO.SaveEdf(image_moving,output_Transform,'moved_rot')
