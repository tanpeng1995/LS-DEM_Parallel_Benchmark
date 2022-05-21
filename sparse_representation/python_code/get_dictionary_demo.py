import glob
import numpy as np
from numpy import savetxt
from numpy import loadtxt
import tifffile as tif
from scipy import signal
import cv2
from scipy.optimize import minimize

img_path = '/Users/tanpeng/Desktop/jupyter_notebooks/'
X = get_compact_X(img_path)
savetxt('../Training/X.csv', X, delimiter=' ')

#D = sparse_coding(X, num_basis = 512, labda = 0.2, num_iters=50, batch_size=500, initD=None)

dic_path = 'Training/D.csv'
D = loadtxt(dic_path, delimiter=' ')

high_dim = D.shape[0]//5
D_high = D[:high_dim,:]
D_low  = D[high_dim:,:]

# normalize the dictionary
# some column is not useful due to zero.
high_norm  = np.sqrt(np.sum(D_high**2, axis=0))
low_norm   = np.sqrt(np.sum(D_low**2, axis=0))
nontrivial = np.intersect1d(np.where(high_norm != 0)[0], np.where(low_norm != 0)[0])

D_high = D_high[:,nontrivial]
D_low  = D_low[:,nontrivial]

D_high = D_high/np.sqrt(np.sum(D_high**2, axis=0))
D_low  = D_low/np.sqrt(np.sum(D_low**2, axis=0))

patch_size = int(np.sqrt(D_high.shape[0]))
labda      = 0.2
upscale    = 2
save_path = 'Training/dictionary'+'_'+str(patch_size)+'_'\
    +str(labda)+'_'+str(upscale)+'_'
savetxt(save_path+'D_high.csv', D_high, delimiter=' ')
savetxt(save_path+'D_low.csv', D_low, delimiter=' ')
