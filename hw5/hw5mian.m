%% hW5

%% Part 1
clear
load('hw5.mat');
N=100;
T=20;
n=20;
beta0=0.1;
sigmabeta=1;
gamma=0;
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
L=sum(logLi,2)

% The calculated log-likelihood function is  -1.5345e+03

%% Part 2
clear
load('hw5.mat');
beta0=0.1;
sigmabeta=1;
gamma=0;
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
L=sum(logLi,2)

% The calculated log-likelihood function is  -1.2385e+03

%% part 3
clear
para=[-0.1; 2; 0];  % initial value
A=[ 0, 0, 0;
    0, -1, 0;
    0, 0, 0];
b= [0; 0; 0];
[para1, min1]= fmincon(@GuaQua, para, A, b)
[para2, min2]= fmincon(@MC, para, A, b)
% Parameter estimates are wrong, wherever I try starting value. but they
% differ with a starting value, which is a bad sign


%  para1 =   0.4772   0.0147   0.0177   min1 = 1.4534e+03
%  para2 =    1.4749   2.5202   0.4761  min2 = 2.6987e+03
% The result is unstable. 


%% Part 4
clear
para= [1; 1; 1; 0.5; 1; 0.5]; % initial value
A= [0,0,0,0,0,0;
    0,0,0,0,0,0;
    0,0,-1,0,0,0;
    0,0,0,0,0,0;
    0,0,0,0,-1,0;
    0,0,0,0,0,0];
b= [0;0;0;0;0;0];
[para3, min3]= fmincon(@MC2, para, A, b)
% these estimates are also wrong, something is wring with your function
% calculation. try to check.

% para3 =
%    1.0001
%    1.0003
%    1.0000
%    0.5002
%    0.9998
%    0.5004
% min3 =
%    1.6338e+03



