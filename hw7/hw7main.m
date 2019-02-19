clear;
tic;
SetupParameters;
% This question is different from the lecture notes. There is no
% competition on quantities. 

V0= repmat(c,1,L);
p0= 10.*ones(L,L);
[Va, pa, itea] =solveMPE(V0,p0);

V0=Va;
p0=pa;
[Va, pa, itea] =solveMPE(V0,p0);


elapsed = toc./60;
disp(sprintf('Elapsed time: %12.4f minutes', elapsed));

figure(1);
mesh(Va);
title('Value Function');

figure(2);
mesh(pa);
title('Policy Function');