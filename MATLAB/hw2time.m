

%Question 4
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
