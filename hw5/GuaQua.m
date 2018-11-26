function L= GuaQua(para)
load('hw5.mat');
beta0=para(1);
sigmabeta=para(2);
gamma=para(3);
N=100;
T=20;
n=20;
[x, w]=qnwnorm(n, beta0, sigmabeta);
F=@(y)((1+exp(-y)).^(-1));
LisWs=zeros(n,N);
for s=1:n
temp=data.X * x(s)+data.Z * gamma ;
Fval=F(temp);
Ls=(Fval.^ data.Y).*(Fval.^(1-data.Y)); 
Lis=prod(Ls,1);
LisWs(s,:)=Lis*w(s);
% logLs= (data.Y).*log(Fval)+ (1-data.Y).*log(1-Fval);
% logLis= sum(logLs, 1);
% LisWs(s,:)= exp(logLis)*w(s);
end
Li= mean(LisWs,1);
logLi= log(Li);
L=-sum(logLi,2);
