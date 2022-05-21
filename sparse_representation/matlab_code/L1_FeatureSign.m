function [x] = L1_FeatureSign(gamma, A, b)
    % this code is reproduced from Yang et al. (2008)

    A           = cast(A, 'double');
    b           = cast(b, 'double');
    eps         = 1e-9;
    x           = zeros(size(A,1), 1);

    grad        = A*sparse(x)+b;
    [max_x, ii] = max(abs(grad).*(x==0));

    while true,
        if grad(ii) > gamma+eps,
            x(ii)   = (gamma-grad(ii))/A(ii, ii);
        elseif grad(ii) < -gamma-eps,
            x(ii)   = (-gamma-grad(ii))/A(ii, ii);
        else
            if all(x==0)
                break;
            end
        end

        while true,
            activated = x~= 0;
            AA = A(activated, activated);
            bb = b(activated);
            xx = x(activated);

            b_new     = -gamma*sign(xx)-bb;
            x_new     = AA\b_new;
            idx       = find(x_new);
            cost_new  = (b_new(idx)/2 + bb(idx))' * x_new(idx) + gamma*sum(abs(x_new(idx)));

            change_signs = find(xx.*x_new <= 0);
            if isempty(change_signs),
                x(activated) = x_new;
                loss         = cost_new;
                break;
            end

            x_min    = x_new;
            cost_min = cost_new;
            d        = x_new-xx;
            t        = d./xx;

            for pos = change_signs',
                x_inter   = xx - d/t(pos);
                x_inter(pos) = 0;
                idx       = find(x_inter);
                cost_temp = (AA(idx, idx)*x_inter(idx)/2 + bb(idx))' * x_inter(idx) + gamma*sum(abs(x_inter(idx)));
                if cost_temp < cost_min,
                    x_min    = x_inter;
                    cost_min = cost_temp;
                end
            end
            x(activated) = x_min;
            loss         = cost_min;
        end

        grad  = A * sparse(x) + b;
        [max_x, ii] = max(abs(grad).*(x==0));
        if max_x <= gamma+eps,
            break;
        end
    end
