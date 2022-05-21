def get_compact_X(img_path, num_basis = 512, labda = 0.15, patch_size = 5, num_patches = 100000, upscale = 2, threshold = 10):
    # randomly sample image patches
    X_high, X_low = rnd_smp_patch(img_path, '*.tif', patch_size, num_patches, upscale)

    # prune patches wiith small variances
    # threshold chosen based on the training data
    X_high, X_low = patch_pruning(X_high, X_low, threshold)

    # train coupled sparse coding dictionary
    # D_high, D_low = train_coupled_dict(X_high, X_low, num_basis, labda, upscale)
    X = train_coupled_dict_before_opt(X_high, X_low, num_basis, labda, upscale)

    # save data
    save_path = 'Training/compact'+'_'+str(num_basis)+'_'\
        +str(labda)+'_'+str(upscale)+'_'
    savetxt(save_path+'X', X, delimiter=' ')

    return X
