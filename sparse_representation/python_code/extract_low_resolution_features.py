def extract_low_resolution_features(low_resolution_img):
    row, col = low_resolution_img.shape
    features = np.zeros((row,col,4))
    # first order gradient filters
    # TODO
    hf1 = np.array([-1,0,1]).reshape(1,-1)
    vf1 = hf1.T
    features[:,:,0] = signal.convolve2d(low_resolution_img, hf1, 'same')
    features[:,:,1] = signal.convolve2d(low_resolution_img, vf1, 'same')

    # second order gradient filters
    # TODO
    hf2 = np.array([1,0,-2,0,1]).reshape(1,-1)
    vf2 = hf2.T
    features[:,:,2] = signal.convolve2d(low_resolution_img, hf2, 'same')
    features[:,:,3] = signal.convolve2d(low_resolution_img, vf2, 'same')
    return features
