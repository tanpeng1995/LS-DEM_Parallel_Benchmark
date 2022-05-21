function [f, g, h] = objective_gradient_hessian(dual_lambda, ZZt, XZt, X, c, trXXt)
    L = size(XZt,1);
    M = length(dual_lambda);
    ZZt_inv = inv(ZZt+diag(dual_lambda));
    if L > M,
        f = -trace(ZZt_inv*(XZt' * XZt)) + trXXt - c*sum(dual_lambda);
    else
        f = -trace(XZt * ZZt_inv * XZt') + trXXt - c*sum(dual_lambda);
    end
    f = -f;

    g = zeros(M,1);
    temp = XZt * ZZt_inv;
    g = sum(temp.^2) - c;
    g = -g;

    h = -2*((temp' * temp).*ZZt_inv);
    h = -h;
