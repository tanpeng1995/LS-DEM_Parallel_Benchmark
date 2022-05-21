def gaussian_kernel(size, sigma=None):
    if sigma == None:
        sigma = size // 6
    center = size // 2
    kernel = np.zeros((size, size))
    for i in np.arange(size):
        for j in np.arange(size):
            diff2 = (i-center)**2 + (j-center)**2
            kernel[i,j] = np.exp(-diff2/sigma**2/2)
    return kernel/np.sum(kernel)

def backprojection(high_resolution_img, low_resolution_img, max_iters):
    row_h, col_h = high_resolution_img.shape
    row_l, col_l = low_resolution_img.shape
    p = gaussian_kernel(5,1)
    p = p**2
    p = p/np.sum(p)

    high_resolution_img = high_resolution_img.astype(float)
    low_resolution_img  = low_resolution_img.astype(float)

    for i in np.arange(max_iters):
        temp = cv2.resize(high_resolution_img, (row_l, col_l), interpolation=cv2.INTER_CUBIC)
        image_diff = low_resolution_img - temp
    
        image_diff = cv2.resize(image_diff, (row_h, col_h), interpolation=cv2.INTER_CUBIC)
        high_resolution_img = high_resolution_img + signal.convolve2d(image_diff, p, 'same')
    return high_resolution_img
