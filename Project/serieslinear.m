%% series estimation linear for muwip0 and muwip1
treated=find(W);
control=find(W0);
sizet=size(treated);
sizec=size(control);
pXhattempt=[ones(sizet(1),1) pXhattemp(treated)];
pXhattempc=[ones(sizec(1),1) pXhattemp(control)];
pXexpand= [ones(N,1) pXhattemp];
pXpXt= pXhattempt'*pXhattempt;
pXpXc= pXhattempc'*pXhattempc;
beta1= pXpXt\pXhattempt'*Y(treated);
beta0= pXpXt\pXhattempc'*Y(control);
muwip1= pXexpand* beta1;
muwip0= pXexpand* beta0;