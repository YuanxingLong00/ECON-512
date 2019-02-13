
clear;
global L luo kapa l c v delta beta lambda CRIT D0 D1 D2; 

L=30;
luo=0.85;
kapa=10;
l=15;
delta=0.03;
beta=1/1.05;
v=10;
lambda=1;
CRIT=1e-10; 

c=zeros(L,1);
c(1:l)=kapa.*([1:l]'.^(log(luo)/ log(2)) );
c(l+1:L)= c(l).*ones(L-l,1);
D0= @(p1,p2)( 1 / (1 + exp(v-p1) +exp(v-p2)) );
D1= @(p1,p2)( exp(v-p1) / (1 + exp(v-p1) +exp(v-p2)) );
D2= @(p1,p2)( exp(v-p2) / (1 + exp(v-p1) +exp(v-p2)) );