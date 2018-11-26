function L= MC(para)
load('hw5.mat');
beta0=para(1);
sigmabeta=para(2);
gamma=para(3);
N=100;
T=20;
S=100; % simulation number
LiNS=zeros(S,N);
F=@(y)((1+exp(-y)).^(-1));
for s=1:S
beta=normrnd(beta0,sigmabeta,1,N);
temp=bsxfun(@times, data.X,beta)+ data.Z * gamma  ;
Fval=F(temp);
Ls=(Fval.^ data.Y).*(Fval.^(1-data.Y)); 
Lis=prod(Ls,1);
LiNS(s,:)=Lis;
end
Li=mean(LiNS,1);
logLi= log(Li);
L=-sum(logLi,2);