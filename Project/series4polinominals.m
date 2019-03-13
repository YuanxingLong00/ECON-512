%% series estimation 4 polinominals for muwip0 and muwip1
treated=find(W);
control=find(W0);
sizet=size(treated);
sizec=size(control);
pXhattempt=[ones(sizet(1),1), pXhattemp(treated), pXhattemp(treated).^2,pXhattemp(treated).^3, pXhattemp(treated).^4];
pXhattempc=[ones(sizec(1),1), pXhattemp(control),pXhattemp(control).^2, pXhattemp(control).^3, pXhattemp(treated).^4];
pXexpand= [ones(N,1), pXhattemp, pXhattemp.^2, pXhattemp.^3, pXhattemp(treated).^4];
pXpXt= pXhattempt'*pXhattempt;
pXpXc= pXhattempc'*pXhattempc;
beta1= pXpXt\pXhattempt'*Y(treated);
beta0= pXpXc\pXhattempc'*Y(control);
muwip1= pXexpand* beta1;
muwip0= pXexpand* beta0;