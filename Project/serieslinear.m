%% series estimation linear for muwip0 and muwip1
pXhreg=[W pXhat];
pXpX=pXhreg'*pXhreg;
beta= pXpX\pXhreg'*Y;
pXhreg1=[ones(N,1) pXhattemp];
pXhreg0=[zeros(N,1) pXhattemp];
muwip1=pXhreg1*beta;
muwip0=pXhreg0*beta;
