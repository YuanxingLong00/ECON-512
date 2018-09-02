function[b_hat, se] = ols(Y, X, hetero)
% Y - dependent variable
% X - regressors (without a constant)
% hetero - 1 if using heteroskedasticity and 0 if using homoskedisticity
    [n, ~] = size(X);

    if numel(Y) ~= n
        error('Incompatible data.') 
    end
    
    % Get estimates
    X1 = [ones(n, 1) X]; XX1 = X1' * X1;
    b_hat = XX1 \ (X1' * Y);
    r_hat = Y - X1 * b_hat;
    
    % Get asyvar.
    switch hetero
        case 1
            X2 = bsxfun(@times, X1, r_hat);
            avar = XX1 \ (X2' * X2) / XX1;
        otherwise
            avar = (r_hat' * r_hat) * inv(XX1) / n; 
    end
 
    se = sqrt(diag(avar));
    
end