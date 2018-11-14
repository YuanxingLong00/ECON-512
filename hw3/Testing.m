%% Testing ?Nelder-Mead simplex method.
clear
t=zeros(1,100);
beta0=[2.5924;-0.0333;0.1162;-0.3572;0.0783;-0.4131];
beta=zeros(100,6);
for i=1:100
load('hw3.mat');
b=beta0+(i-1)*[0.01;0.01;0.01;0.01;0.01;0.01];
tic;
beta(i,:) =fminsearch(@ml,b);
t(i)=toc;
end
time=mean(t);
time=mean(t(1:32));


%% Testing Broyden's method. 
clear
beta0= [2.5924;-0.0333;0.1162;-0.3572;0.0783;-0.4131];
t=zeros(1,100);
beta_rec=zeros(100,6);

for i=1:100
    beta=beta0+(i-1)*[0.01;0.01;0.01;0.01;0.01;0.01];
    tic;
fVal = FocML(beta);
iJac = inv(JacFoc(beta));
maxit = 1000;
tol = 1e-6;
for iter = 1:maxit
    fnorm = norm(fVal);
    if norm(fVal) < tol
        beta_rec(i,:)=beta;
        fprintf('testing %d  Convergence\n', i);
        break
    end
    d = - (iJac * fVal);
    beta = beta +d;
    fOld = fVal;
    fVal = FocML(beta);
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end
t(i)=toc;
end
time=mean(t);
time=mean(t(1:4));

%% Testing Nonlinear LS lsqnonlin with option trust-region-reflective algorithm.
clear 
t=zeros(1,100);
beta0=[2.5924;-0.0333;0.1162;-0.3572;0.0783;-0.4131];
beta_rec=zeros(100,6);
for i=1:100
load('hw3.mat');
b=beta0+(i-1)*[0.01;0.01;0.01;0.01;0.01;0.01];
tic;
fun=@(beta) y- exp(X*[beta(1);beta(2);beta(3);beta(4);beta(5);beta(6)]);
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
tic;
beta_rec(i,:)=lsqnonlin(fun,b,[],[],options);
t(i)=toc;
end
time=mean(t);
time=mean(t(1:99))

