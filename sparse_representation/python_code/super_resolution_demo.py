# Reference:
# J. Yang et al. Image super-resolution via sparse representation. IEEE
# Transactions on Image Processing, Vol 19, Issue 11, pp2861-2873, 2010
#
# This code is re-written by Peng TAN (Berkeley)
# source code is referred from Prof. JianChao Yang, ECE department, UIUC
import glob
import numpy as np
from numpy import savetxt
from numpy import loadtxt
import tifffile as tif
from scipy import signal
import cv2

labda       = 0.2   # sparsity regularization
overlap     = 4     # the more overlap the better (patch size 5x5)
upscale     = 2     # scaling factor, depending on the trained dictionary
max_iters   = 20    # if 0, do not use backprojection
patch_size  = 5

D_high_path = 'Training/dictionary'+'_'+str(patch_size)+'_'\
    +str(labda)+'_'+str(upscale)+'_D_high.csv'
D_low_path  = 'Training/dictionary'+'_'+str(patch_size)+'_'\
    +str(labda)+'_'+str(upscale)+'_D_low.csv'
D_high = loadtxt(D_high_path, delimiter=' ')
D_low  = loadtxt(D_low_path, delimiter=' ')

low_resolution_img_path = 'low_resolution_img.tif'
low_resolution_img = tif.imread(low_resolution_img_path)

# image super-resolution based on sparse representation
# cross validation TODO
high_resolution_img = img_super_resolution(low_resolution_img, upscale, D_high, D_low, labda, overlap)
high_resolution_img = backprojection(high_resolution_img, low_resolution_img, max_iters)
