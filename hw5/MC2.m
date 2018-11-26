function L= MC2(para)
load('hw5.mat');
meanmu=[para(1); para(2)];
Sigma= [para(3), para(4);
        para(4),para(5)];
gamma=para(6);

% This step is to guarantee Sigma is positive semidefinite.
if para(3)*para(5)-para(4)^2 <0 || para(3)+ para(5) < 0
    L=-Inf;
else
    
    
N=100;
T=20;
S=100; % simulation number
LiNS=zeros(S,N);
one= ones(T,N);
F=@(y)((1+exp(-y)).^(-1));
for s=1:S
mu= mvnrnd(meanmu, Sigma,N);
beta=mu(:,1)';
u=mu(:,2)';
% beta=normrnd(beta0,sigmabeta,1,N);
temp=bsxfun(@times, data.X,beta)+ data.Z * gamma + bsxfun(@times, one, u);
Fval=F(temp);
Ls=(Fval.^ data.Y).*(Fval.^(1-data.Y)); 
Lis=prod(Ls,1);
LiNS(s,:)=Lis;
end
Li=mean(LiNS,1);
logLi= log(Li);
L=-sum(logLi,2);

end