function [Z] = L1_FeatureSign_Setup(X, D, labda)
    [M, N] = size(X);
    K      = size(D,2);
    Z      = sparse(K, N);
    A      = D' * D;
    for i  = 1:N,
        if mod(i,1000) == 0, disp(i); end
        b  = -D' * X(:,i);
        Z(:,i) = L1_FeatureSign(labda, A, b);
    end
