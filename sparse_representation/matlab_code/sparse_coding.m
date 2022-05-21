function [D, Z] = sparse_coding(X, num_basis, labda, num_iters, batch_size)
    [M, N] = size(X);
    K      = num_basis;
    X      = cast(X, 'double');
    if isempty(batch_size) batch_size = N; end
    D      = rand(M, K) - 0.5;
    D      = D - repmat(mean(D,1), M, 1);
    D      = D ./ repmat(sqrt(sum(D.^2)), M, 1);

    count = 0;
    % optimization loop
    while count < num_iters
        %disp(count);
        count = count + 1;
        idx   = randperm(N);
        for i = 1:(N/batch_size)
            batch_idx = idx((1:batch_size)+batch_size*(i-1));
            X_batch   = X(:,idx);
            Z         = L1_FeatureSign_Setup(X_batch, D, labda);
            D         = L2_Lagrange_Dual(X_batch, Z, 1);
        end
    end

    writematrix(D, 'Training/D.csv', 'Delimiter', 'space')
