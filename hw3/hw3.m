% all is great

%% HW3 
%% Q1 MLE using Nelder_Mead simplex 
clear
load('hw3.mat');
beta0=[1;1;1;1;1;1];
options = optimset('PlotFcns',@optimplotfval, 'Display','iter');
beta =fminsearch(@ml,beta0, options);
beta

%% Q2 Use Broyden's method to solve ML estimator
clear
beta= [2.5924;-0.0333;0.1162;-0.3572;0.0783;-0.4131];
fVal = FocML(beta);
iJac = inv(JacFoc(beta));

maxit = 100;
tol = 1e-6;
for iter = 1:maxit
    fnorm = norm(fVal);
    g=sprintf('%.4f', beta);
    fprintf('iter %d: beta= %s, norm(f(x)) = %.8f\n', iter, g, norm(fVal));
    if norm(fVal) < tol
        break
    end
    d = - (iJac * fVal);
    beta = beta +d;
    fOld = fVal;
    fVal = FocML(beta);
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end

%% Q3 Use Nonlinear LS to estimate beta
clear 
load('hw3.mat');
fun=@(beta) y- exp(X*[beta(1);beta(2);beta(3);beta(4);beta(5);beta(6)]);
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
beta0=[2.5; -0.0; 0.1; -0.4; 0.1; -0.5];
beta=lsqnonlin(fun,beta0,[],[],options)

%% Q4 ?Estimate the parameter vector $\beta$ using the nonlinear least squares estimator computed using the Nelder-Mead simplex method. 
clear
beta0=[2.5122; -0.0384; 0.1141; -0.2799;  0.0677; -0.3697];
options = optimset('PlotFcns',@optimplotfval, 'Display','iter');
beta =fminsearch(@ml,beta0, options);
beta

% Where is question 5? 
% also, all the answers don't look right. did you try to squeeze the
% stopping criteria?