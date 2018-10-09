function Jac=JacFoc(beta)
load('hw3.mat');
temp= -exp(X*[beta(1);beta(2);beta(3);beta(4);beta(5);beta(6)]);
temp1=(temp.*X)';
Jac=temp1*X;

