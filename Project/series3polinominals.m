%% series estimation 3 polinominals for muwip0 and muwip1
pXhreg=[W pXhat pXhat.^2 pXhat.^3];
pXpX=pXhreg'*pXhreg;
beta= pXpX\pXhreg'*Y;
pXhreg1=[ones(N,1) pXhattemp.^2 pXhattemp.^3];
pXhreg0=[zeros(N,1) pXhattemp.^2 pXhattemp.^3];
muwip1=pXhreg1*beta;
muwip0=pXhreg0*beta;