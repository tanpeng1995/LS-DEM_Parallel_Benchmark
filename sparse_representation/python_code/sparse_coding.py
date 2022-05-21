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



def L1_FeatureSign_Setup(X, D, labda):
    """
    minimize ||X[:,j] - D*Z[:,j]||_2^2 + 2*sigma^2*beta |Z[:,j]|_1
    2*sigma^2*beta = labda
    since z = Z[:,j] could be separately optimized
    """
    M, N = X.shape
    K = D.shape[1]
    Z = np.zeros((K,N))
    A = D.T @ D
    for i in np.arange(N):
        b = -D.T @ X[:,i]
        #sklearn
        Z[:,i] = L1_FeatureSign(labda, A, b)
    return Z

def object_func(dual_lambda, ZZt, XZt, X, c, trXXt):
    # objective function at given dual_lambda
    L = XZt.shape[0]
    M = len(dual_lambda)
    # TODO
    ZZt_inv = np.linalg.inv(ZZt+np.diag(dual_lambda))
    if L > M:
        f = -np.trace(ZZt_inv @ (XZt.T @ XZt)) + trXXt - c*np.sum(dual_lambda)
    else:
        f = -np.trace(XZt @ ZZt_inv @ XZt.T) + trXXt - c*np.sum(dual_lambda)
    return -f

def gradient_func(dual_lambda, ZZt, XZt, X, c, trXXt):
    L = XZt.shape[0]
    M = len(dual_lambda)
    # TODO
    ZZt_inv = np.linalg.inv(ZZt+np.diag(dual_lambda))
    # gradient of the function at given dual_lambda
    g = np.zeros((M,1))
    temp = XZt @ ZZt_inv
    g = np.sum(temp**2, axis=0)-c
    return -g

def hessian_func(dual_lambda, ZZt, XZt, X, c, trXXt):
    L = XZt.shape[0]
    M = len(dual_lambda)
    ZZt_inv = np.linalg.inv(ZZt+np.diag(dual_lambda))
    # Hessian evaluated at given dual_lambda
    temp = XZt @ ZZt_inv
    h = -2 * (temp.T @ temp) * ZZt_inv
    return -h

def L2_Lagrange_Dual(X, Z, c):
    M, N = X.shape
    K = Z.shape[0]

    ZZt = Z @ Z.T
    XZt = X @ Z.T

    # arbitrary initialization dual_lambda as Gaussian
    dual_lambda = 10*np.abs(np.random.rand(K))
    trXXt = np.sum(X**2)
    options = {'disp':True, 'maxiter': 100}
    # those three are good, if warning says Desired error not necessarily achieved due to precision loss.
    # try another method
    res = minimize(object_func, dual_lambda, args=(ZZt, XZt, X, c, trXXt), jac = gradient_func, method='CG', options=options)
    #res = minimize(object_func, dual_lambda, args=(ZZt, XZt, X, c, trXXt), jac = gradient_func, hess = hessian_func, method='trust-ncg', options=options)
    #res = minimize(object_func, dual_lambda, args=(ZZt, XZt, X, c, trXXt), jac = gradient_func, hess = hessian_func, method='Newton-CG', options=options)
    dual_lambda = res.x
    # TODO
    Dt = np.linalg.inv(ZZt+np.diag(dual_lambda)) @ XZt.T
    D  = Dt.T
    return D

def sparse_coding(X, num_basis, labda, num_iters=50, batch_size=500, initD=None):
    """
    Regularized Sparse Coding
    X:          preprocessed np.ndarray, dimension: MxN
    num_basis:  number of basis. D.shape[1]
    gamma:      sparsity regularization
    num_iters:  number of iterations
    batch_size: batch size
    initD:      initial dictionary

    D:          learned dictionary, dimension: MxK, where K << N, K=num_basis
    Z:          sparse code, dimension: KxN

    this function solve:
        minimize_{D,Z} ||X-DZ||_2^2 + lambda*|Z|_1
        s.t. column norm of D <= c
    This is not convex in both D and Z, but is convex in one of them with the other fixed.
    """
    M, N = X.shape # M: patch_size, N: num_patches
    K    = num_basis
    X = X.astype(float)
    if batch_size == None: batch_size = N
    if initD == None:
        # Initialize D with a Gaussian random matrix
        D = np.random.rand(M, K)-0.5
        D = D - np.mean(D, axis=0)
        # each column is unit normalized
        D = D/np.sqrt(np.sum(D**2, axis=0))
    else:
        D = initD

    # optimization loop:
    for iter in np.arange(num_iters):
        print('{}-th iterations for sparse coding'.format(iter))
        idx = np.arange(N)
        np.random.shuffle(idx)

        for i in np.arange(N//batch_size):
            batch_idx = idx[np.arange(batch_size)+batch_size*(i-1)]
            X_batch   = X[:,batch_idx]

            # fix D, update Z
            # Z = argmin_Z ||X-DZ||_2^2 + labda*|Z|_1
            # for the paper referenced, labda = sigma^2 * beta
            # this is provided by Honglak Lee. et al (2017), Efficient sparse coding algorithms
            Z = L1_FeatureSign_Setup(X_batch, D, labda)

            # fix Z, update D
            # c = 1 sum_j ||D_i^j|| \leq 1, i = 1,2,...,k, column norm constraints
            D = L2_Lagrange_Dual(X_batch, Z, 1)

    return D, Z
