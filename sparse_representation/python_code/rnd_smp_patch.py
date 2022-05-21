def normalize_image(img, L):
    _min, _max = img.min(), img.max()
    img = ((img-_min)/(_max-_min) * (L-1)).astype(np.uint8)
    return img

def sample_patches(img, patch_size, num_copies, upscale):
    high_resolution_img = img
    high_resolution_img = high_resolution_img.astype(float)
    m, n = high_resolution_img.shape
    # generate low resolution counterparts
    low_resolution_img  = cv2.resize(\
        high_resolution_img, (m//upscale, n//upscale),interpolation=cv2.INTER_CUBIC)
    low_resolution_img  = cv2.resize(\
        low_resolution_img, (m, n), interpolation=cv2.INTER_CUBIC)
    low_resolution_img = low_resolution_img.astype(float)

    x = np.arange(m-2*patch_size)+patch_size
    y = np.arange(n-2*patch_size)+patch_size
    np.random.shuffle(x)
    np.random.shuffle(y)

    X, Y = np.meshgrid(x,y)
    xrow, ycol = X.reshape(-1), Y.reshape(-1)
    if num_copies < len(xrow):
        xrow = xrow[:num_copies]
        ycol = ycol[:num_copies]
    else:
        num_copies = len(xrow)

    # initialize output
    x_high = np.zeros((patch_size**2, num_copies))
    x_low  = np.zeros((4*patch_size**2, num_copies))

    # compute the first and second order gradients
    hf1 = np.array([-1,0,1]).reshape(1,-1)
    vf1 = hf1.T
    hf2 = np.array([1,0,-2,0,1]).reshape(1,-1)
    vf2 = hf2.T

    # get low_resolution_img features
    low_resolution_feature1 = signal.convolve2d(low_resolution_img, hf1, 'same')
    low_resolution_feature2 = signal.convolve2d(low_resolution_img, vf1, 'same')
    low_resolution_feature3 = signal.convolve2d(low_resolution_img, hf2, 'same')
    low_resolution_feature4 = signal.convolve2d(low_resolution_img, vf2, 'same')

    # collect patches from sample
    for i in np.arange(num_copies):
        row, col = xrow[i], ycol[i]
        Hpatch  = high_resolution_img[row:row+patch_size, col:col+patch_size].reshape(-1)
        Lpatch1 = low_resolution_feature1[row:row+patch_size, col:col+patch_size].reshape(-1)
        Lpatch2 = low_resolution_feature2[row:row+patch_size, col:col+patch_size].reshape(-1)
        Lpatch3 = low_resolution_feature3[row:row+patch_size, col:col+patch_size].reshape(-1)
        Lpatch4 = low_resolution_feature4[row:row+patch_size, col:col+patch_size].reshape(-1)
        Lpatch  = np.concatenate([Lpatch1,Lpatch2,Lpatch3,Lpatch4],axis=0)
        x_high[:,i] = Hpatch-np.mean(Hpatch)
        x_low[:,i]  = Lpatch

    return x_high, x_low

def rnd_smp_patch(img_path, type, patch_size, num_patches, upscale):
    # get all training images name
    img_list = glob.glob(img_path+type) # type = '*.tif'
    # get total number of images being considered
    img_num = len(img_list)
    # initialize number of copies for each image
    # depends on its size
    num_copies_img = np.zeros(img_num)

    # read images and determine number of copies for each image
    # this number is proportional to total number of patches
    for i in np.arange(img_num):
        img = tif.imread(img_list[i])
        num_copies_img[i] = np.prod(img.shape)
    num_copies_img = np.floor(num_copies_img*num_patches/np.sum(num_copies_img)).astype(np.int)

    # initialize output
    X_high = []
    X_low  = []

    for i in np.arange(img_num):
        num_copies = num_copies_img[i]
        img = tif.imread(img_list[i])
        img = normalize_image(img, 256)
        x_high, x_low = sample_patches(img, patch_size, num_copies, upscale)
        X_high.append(x_high)
        X_low.append(x_low)

    # assemble a numpy ndarray
    X_high = np.concatenate(X_high, axis=1)
    X_low  = np.concatenate(X_low, axis=1)

    # save data
    save_path = 'Training/rnd_patches'+str(patch_size)+'_'+str(upscale)+'_'
    savetxt(save_path+'X_high.csv', X_high, delimiter=' ')
    savetxt(save_path+'X_low.csv', X_low, delimiter=' ')

    return X_high, X_low

def patch_pruning(X_high, X_low, threshold):
    vars = np.var(X_high, axis=0)
    idx  = vars > threshold
    X_high = X_high[:,idx]
    X_low  = X_low[:,idx]
    return X_high, X_low
