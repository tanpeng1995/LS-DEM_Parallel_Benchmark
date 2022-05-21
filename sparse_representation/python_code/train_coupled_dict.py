def train_coupled_dict(X_high, X_low, num_basis, labda, upscale):
    """
    solve:
    D = argmin ||X-DZ||_2^2 + labda * |Z|_1
    s.t ||D_i||_2^2 <= 1, for i = 1, 2, 3, ..., num_basis
    """
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

    # dictionary training
    D, Z = sparse_coding(X, num_basis, labda) # X = DZ, num_basis = D.shape[1]
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

    return D_high, D_low

def train_coupled_dict_before_opt(X_high, X_low, num_basis, labda, upscale):
    """
    Unfortunately, I did not find python equivalent to solve this function
    I have to go to MATLAB and carry the result back.
    solve:
    D = argmin ||X-DZ||_2^2 + labda * |Z|_1
    s.t ||D_i||_2^2 <= 1, for i = 1, 2, 3, ..., num_basis
    """
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

    return X


def train_coupled_dict_after_opt(D, Z):
    # dictionary training
    D, Z = sparse_coding(X, num_basis, labda) # X = DZ, num_basis = D.shape[1]
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

    return D_high, D_low
