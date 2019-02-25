function [V,p,iter] = solveMPE(V0, p0)
% using global variables in matlab is slowing it down a lot. 
lambda=1;
CRIT=1e-5; 


V=V0;
p=p0;
iter=0;

for i=1:1500
    [W0,W1,W2]=GetW(V);
    np=GetP(W0,W1,W2,p); 
    nV=GetV(p,np,W0,W1,W2);
    check = max( max(max(abs((nV-V)./(1+nV)))),  max(max(abs((np-p)./(1+nV))))     )
    V = lambda.* nV + (1-lambda).*V;
    p = lambda.* np + (1-lambda).*p;
   if check<CRIT
       fprintf('Convergence Achieved');
        break;
   end
    iter=i
end

 end