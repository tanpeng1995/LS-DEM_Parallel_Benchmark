import numpy as np
def normalize_image(img, L):
    _min, _max = img.min(), img.max()
    img = ((img-_min)/(_max-_min) * (L-1)).astype(np.uint8)
    return img

def patch_pruning(X_high, X_low, threshold):
    """
    Some patchs trivial and do not contain too much information, delete them according to threshold.
    There is no code work for this part.
    This function provides data pre-processing tools.
    """
    vars = np.var(X_high, axis=0)
    idx  = vars > threshold
    X_high = X_high[:,idx]
    X_low  = X_low[:,idx]

    high_dim = X_high.shape[0]
    low_dim  = X_low.shape[0]

    # should pre-normalize X_high and X_low
    high_norm  = np.sqrt(np.sum(X_high**2, axis=0))
    low_norm   = np.sqrt(np.sum(X_low**2, axis=0))
    nontrivial = np.intersect1d(np.where(high_norm != 0)[0], np.where(low_norm != 0)[0])

    X_high = X_high[:,nontrivial]
    X_low  = X_low[:,nontrivial]

    X_high = X_high/np.sqrt(np.sum(X_high**2, axis=0))
    X_low  = X_low/np.sqrt(np.sum(X_low**2, axis=0))

    # joint learning of the dictionary
    X = np.concatenate([np.sqrt(high_dim)*X_high, np.sqrt(low_dim)*X_low], axis=0)
    X_norm = np.sqrt(np.sum(X**2, axis=0))
    X = X[:,X_norm > 1e-5]
    X = X/np.sqrt(np.sum(X**2, axis=0))
    
    # save data
    # save_path = '../Training/rnd_patches'+str(patch_size)+'_'+str(upscale)+'_'
    # savetxt(save_path+'X.csv', X, delimiter=' ')

    return X

def L1_FeatureSign(gamma, A, b):
    """
    The detail of the algorithm is described in the following paper:
    'Efficient Sparse Coding Algorithms', Honglak Lee, Alexis Battle, Rajat Raina, Andrew Y. Ng,
    Advances in Neural Information Processing Systems (NIPS) 19, 2007

    minimize 0.5*x.T @ A @ x + b.T @ x + gamma * |x|
    """
    A = A.astype(float)
    b = b.astype(float)
    x = np.zeros(A.shape[0])
    eps = 1e-9
    grad = A @ x + b
    ii = np.argmax(np.abs(grad*(x==0)))

    while True:
        if grad[ii] > gamma+eps:
            x[ii] = (gamma-grad[ii])/A[ii,ii]
        elif grad[ii] < -gamma-eps:
            x[ii] = (-gamma-grad[ii])/A[ii,ii]
        else:
            if np.all(x==0):
                break

        while True:
            # consider active set
            activated = x != 0
            AA = A[activated,:][:,activated]
            bb = b[activated]
            xx = x[activated]
            # new b based on unchanged sign
            # Ax + b + gamma * sign(x) = 0
            b_new = -gamma*np.sign(xx)-bb
            # analytical solution
            x_new = np.linalg.inv(AA)@b_new
            idx   = x_new != 0
            cost_new  = (b_new[idx]/2 + bb[idx]).T @ x_new[idx] + \
                                gamma*np.sum(np.abs(x_new[idx]))
            change_signs = np.where(xx*x_new <= 0)[0]
            # if no sign change, x_new is optimum since it's analytical solution
            if len(change_signs) == 0:
                x[activated] = x_new
                loss         = cost_new
                break
            # find the best interpolation solution x_inter between x_new and xx
            # x_inter is improved compared with xx, because of convexity
            x_min    = x_new
            cost_min = cost_new
            d        = x_new-xx
            t        = d/xx
            for pos in change_signs:
                x_inter = xx - d/t[pos] # interpolating at pos-th point
                x_inter[pos] = 0        # make sure it is zero
                idx = x_inter != 0
                # cost of x_inter
                cost_temp = (AA[idx,:][:,idx]@x_inter[idx]/2+bb[idx]).T@x_inter[idx] + \
                        gamma*np.sum(np.abs(x_inter[idx]))
                if cost_temp < cost_min:
                    x_min = x_inter
                    cost_min = cost_temp
            # update x and loss
            x[activated] = x_min
            loss         = cost_min

        grad  = A @ x + b
        ii    = np.argmax(np.abs(grad*(x==0)))
        max_x = np.abs(grad[ii])
        if max_x <= gamma+eps:
            break

    return x
