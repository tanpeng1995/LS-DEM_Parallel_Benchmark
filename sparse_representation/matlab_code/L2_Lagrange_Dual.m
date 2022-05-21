function [D] = L2_Lagrange_Dual(X, Z, c)
    [M, N]      = size(X);
    K           = size(Z,1);

    ZZt         = Z*Z';
    XZt         = X*Z';
    dual_lambda = 10*abs(rand(K,1));
    trXXt       = sum(sum(X.^2));
    lb          = zeros(size(dual_lambda));
    options     = optimset('GradObj','on', 'HessFcn','objective', 'Display', 'iter', 'TolFun', 0.0001,'TolX', 0.0001);
    x           = fmincon(@(x) objective_gradient_hessian(x, ZZt, XZt, X, c, trXXt), dual_lambda, [], [], [], [], lb, [], [], options);
    dual_lambda = x;
    Dt          = (ZZt+diag(dual_lambda)) \ XZt';
    D           = Dt';
