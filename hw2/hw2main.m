%% Econ 512 HW2
% Question 2 use Broyden's method to solve for Nash pricing Eqm
% find initial values for jacobian
clear
syms x y
g=@(x,y)[1-x+x*exp(2-x)/(1+exp(2-x)+exp(2-y)),1-y+y*exp(2-y)/(1+exp(2-x)+exp(2-y))];
jacob=jacobian(g,[x,y])
% From above calculations, I find the Jacobian matrix and use jacob
% function to evaluate it at (x,y). 
% Next specify the initial values for p, jacobian and fVal
clear
f=@(x)[1-x(1)+x(1)*exp(2-x(1))/(1+exp(2-x(1))+exp(2-x(2))),1-x(2)+x(2)*exp(2-x(2))/(1+exp(2-x(1))+exp(2-x(2)))];
p=[1.5;1.5] % why these parameter values? different people use these parameters. 
fVal=f(p)'
Jac=jacob(p)  % jacob function is defined seperately by using the result above
iJac=inv(Jac)
% Now do the Broyden's iterations:
maxit=100;
tol=1e-6;
for iter=1:maxit
    fnorm = norm(fVal);
    fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n', iter, p(1), p(2), norm(fVal));
    if fnorm < tol
        break
    end
    d = - (iJac * fVal);
    p=p+d;
    fOld=fVal;
    fVal=f(p)';
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end

%% Question 3 
clear
p=[1.5;1.5]; % initial value for p
pA=p(1);
pB=p(2);
pAOld=1.4;
pBOld=1.4;
hOld=h(pAOld, pB);
%secant iterations
tol = 1e-8;
maxit = 100;
for iter=1:maxit    
for i =1:maxit
    hVal = h(pA,pB);
    fprintf('iter %d: pA = %.8f, pB= %.8f, h(pA, pB) = %.8f\n', i, pA, pB, hVal);
    if abs(hVal) < tol
        break
    else
        pANew = pA - ( (pA - pAOld) / (hVal - hOld) )* hVal;
        pAOld = pA;
        pA = pANew;
        hOld = hVal;       
    end
end
pBlag=pB;
for j =1:maxit
    hVal = h(pB,pA);% use the symmetry of FOC for firm A and B.
    fprintf('iter %d: pB = %.8f, pA= %.8f, h(pB, pA) = %.8f\n', j, pB, pA, hVal);
    if abs(hVal) < tol
        break
    else
        pBNew = pB - ( (pB - pBOld) / (hVal - hOld) )* hVal;
        pBOld = pB;
        pB = pBNew;
        hOld = hVal;       
    end
end
if abs(pB-pBlag)<tol
    fprintf('iter %d: pA = %.8f, pB= %.8f\n', iter, pA, pB);
    break
end
end

%% Question 4
clear 
% initial guess for pA and pB
pA=1.5;
pB=1.5;
tol = 1e-8;
maxit = 100;
for iter=1:maxit  
    pAlag=pA;
    pBlag=pB;
    [pA,pB]=update(pA,pB);
    if abs(pA-pAlag)<tol && abs(pB-pBlag)<tol
        fprintf('iter %d: pA = %.8f, pB= %.8f\n', iter, pA, pB);
        break
    end
end


%% Question 5
% need to modify it.
clear
pA=0:.2:3;
pB=0:.2:3;
VB=0:.2:3;
maxit=100;
tol=1e-6;
for i=0:15
 vB=VB(i+1);
f=@(x)[1-x(1)+x(1)*exp(2-x(1))/(1+exp(2-x(1))+exp(vB-x(2))),1-x(2)+x(2)*exp(vB-x(2))/(1+exp(2-x(1))+exp(vB-x(2)))];
p=[1.5;1.5];
fVal=f(p)';
Jac=jacob(p);  % jacob function is defined seperately by using the result above
iJac=inv(Jac);
% Now do the Broyden's iterations:
for iter=1:maxit
    fnorm = norm(fVal);
    if fnorm < tol
        fprintf('iter %d: pA = %f, pB = %f, norm(f(x)) = %.8f\n', iter, p(1), p(2), norm(fVal));
        break
    end
    d = - (iJac * fVal);
    p=p+d;
    fOld=fVal;
    fVal=f(p)';
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end
pA(i+1)=p(1);
pB(i+1)=p(2);
end
plot(VB,pA,'*',VB,pB,'*')

