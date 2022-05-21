def lin_scale(hPatch, mNorm):
    hNorm = np.linalg.norm(hPatch)
    if hNorm > 0:
        scale = 1.2 * mNorm / hNorm
        hPatch = scale * hPatch
    return hPatch

def normalize_image(img, L):
    _min, _max = img.min(), img.max()
    img = ((img-_min)/(_max-_min) * (L-1)).astype(np.uint8)
    return img

def img_super_resolution(low_resolution_img, upscale, D_high, D_low, labda, overlap):
    # normalize the dictionary
    D_high = D_high/np.sqrt(np.sum(D_high**2, axis=0))
    D_low  = D_low/np.sqrt(np.sum(D_low**2, axis=0))

    # get patch_size
    patch_size = int(np.sqrt(D_high.shape[0]))

    # bicubic interpolation of the low_resolution_img
    m, n = low_resolution_img.shape
    # TODO
    medium_resolution_img = cv2.resize(low_resolution_img, \
        (m*upscale, n*upscale), interpolation=cv2.INTER_CUBIC)
    M, N = medium_resolution_img.shape

    # initialize high_resolution_img
    high_resolution_img = np.zeros(medium_resolution_img.shape)
    count_map = np.zeros(medium_resolution_img.shape)

    # extract low_resolution_img features
    features = extract_low_resolution_features(medium_resolution_img)

    # patch index for sparse recovery, avoid boundary
    p = patch_size//2
    gridx = np.array(list(range(p,M-patch_size-p,patch_size-overlap))+[M-patch_size-p])
    gridy = np.array(list(range(p,N-patch_size-p,patch_size-overlap))+[N-patch_size-p])

    A = D_low.T @ D_low

    # loop to recover each high resolution patch
    for i in np.arange(len(gridx)):
        for j in np.arange(len(gridx)):
            patch_idx = i*len(gridy)+j
            xx = gridx[i]
            yy = gridy[j]

            # column feature
            mPatch = medium_resolution_img[xx:xx+patch_size, yy:yy+patch_size].reshape(-1)
            mMean  = np.mean(mPatch)
            mPatch = mPatch - mMean
            mNorm  = np.linalg.norm(mPatch)

            mPatchFea = features[xx:xx+patch_size, yy:yy+patch_size, :].reshape(-1)
            mPatchFea = mPatchFea - np.mean(mPatchFea)
            mFeaNorm  = np.linalg.norm(mPatchFea)

            y = mPatchFea / mNorm if mFeaNorm > 1 else mPatchFea
            b = -D_low.T @ y

            # solve for sparce coefficient using feature sign
            w = L1_FeatureSign(labda, A, b)

            # recover high resolution patch and scale the contrast
            hPatch = D_high @ w
            hPatch = lin_scale(hPatch, mNorm)
            hPatch = hPatch.reshape(patch_size, patch_size)
            hPatch = hPatch + mMean

            high_resolution_img[xx:xx+patch_size, yy:yy+patch_size] += hPatch
            count_map[xx:xx+patch_size, yy:yy+patch_size] += 1

    # fill in the empty with bicubic interpolation
    high_resolution_img[count_map < 1] = medium_resolution_img[count_map < 1]
    count_map[count_map < 1] = 1
    high_resolution_img = high_resolution_img / count_map
    high_resolution_img = normalize_image(high_resolution_img, 256)

    return high_resolution_img
