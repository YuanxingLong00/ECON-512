clear;
tic;
SetupParameters; 
% This question is different from the lecture notes. There is no
% competition on quantities. 
L=30;

V0= 20.*ones(L,L);
p0= 10.*ones(L,L);
[Va, pa, itea] =solveMPE(V0,p0);

elapsed = toc./60;
fprintf('Elapsed time: %12.4f minutes', elapsed);

figure(1);
mesh(Va);
title('Value Function');

figure(2);
mesh(pa);
title('Policy Function');

%% plot the state distribution after 10, 20 and 30 periods
plot1=zeros(L,L);
plot2=zeros(L,L);
plot3=zeros(L,L);
w0=[1 1];
for i=1:1000
W10=getst(w0,10,pa);
plot1(W10)=plot1(W10)+1;

W20=getst(w0,20,pa);
plot2(W20)=plot2(W20)+1;

W30=getst(w0,30,pa);
plot3(W30)=plot3(W30)+1;
end
figure(3);
mesh(plot1);
title('State Distribution after 10 periods');

figure(4);
mesh(plot2);
title('State Distribution after 20 periods');

figure(5);
mesh(plot3);
title('Sate Distribution after 30 periods');

%% get the stationary distribution of states
plot4=zeros(L,L);
w0=[15,15];
for i=1:1000
    W10=getst(w0,10000,pa);
plot4(W10)=plot4(W10)+1;
end
figure(6);
mesh(plot4);
title('Stationary Sate Distribution');

%% The distribution looks weired. 
%% Then best way is to calculate the transition probability matrix and then 
%% matrix multiplication to calculate the future distribution. 
%% The stationary distribution is just the solution to y=Py where P is transition matrix. 
