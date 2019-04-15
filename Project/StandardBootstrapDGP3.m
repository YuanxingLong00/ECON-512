%%%% Standard Bootstrap method for propensity score matching estimator.
%%%% DGP3
%%%% H0: tau<=5 vs H1: tau>5
% clear;
rng(246);
% N=2000;
alpha=0.05;
B=500;
R=200;
Rej=0;
M=0;

for r=1:R
% generate data
N1=M;
N0=M;
while N1<= M || N0<= M
X1= rand(N,1)-0.5*ones(N,1);
X2= rand(N,1)-0.5*ones(N,1);
X=[X1 X2];
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= 3*X1-3*X2+U0;
Y1= 5+5*X1+X2+U1;
tau=5;
pX= (exp(X1+7*X2))./(ones(N,1)+exp(X1+7*X2));
W= binornd(1,pX);
W0=ones(N,1)-W;
N1=sum(W);
N0=N-N1;
if  N1<= M || N0<=M
     continue
 end
 break
end


Y= W.*Y1+(ones(N,1)-W).*Y0;
% Estimate theta
F= @(theta) (exp(theta(1)*X1+theta(2)*X2)./(ones(N,1)+exp(theta(1)*X1+theta(2)*X2)));
L= @(theta) (-sum(W.*log(F(theta))+(ones(N,1)-W).*log(ones(N,1)-F(theta))));
theta0=[1;7];
thetahat= fminsearch(L,theta0);
pXhat=F(thetahat); % Estimated Propensity Score

% first, construct estimate for ATE tau based on estimated propensity score
Ybar1=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
for i=1:N
    if W(i)==1
        DPS= abs( pXhat(i)*W0-pXhat.*W0 );
        DPS(W==1)=2;
        [mi,ind]=min(DPS); 
        Ybar1(i)=Y(ind);
    else 
        DPS= abs(pXhat(i)*W-pXhat.*W);
        DPS(W==0)=2;
        [mi,ind]=min(DPS);
        Ybar1(i)=Y(ind);
    end
end
tauhat= sum( ( 2*W-ones(N,1) ).*(Y-Ybar1) )/N; % tauhat based on estimated propensity score 

Tdistr=zeros(B,1);
for b=1:B
N1s= M;
N0s= M;
while N1s<= M || N0s<=M % Eliminate unbalanced bootstrap sample 
  S= unidrnd(N,N,1);
    Xs= X(S,:);
    Ys= Y(S);
    Ws= W(S);
    W0s= ones(N,1)-Ws;
 N1s= sum(Ws);
 N0s= N- N1s;
 if  N1s<= M || N0s<=M
     continue
 end
 break
end
X1s= Xs(:,1);
X2s= Xs(:,2);

F= @(theta) (exp(theta(1)*X1s+theta(2)*X2s)./(ones(N,1)+exp(theta(1)*X1s+theta(2)*X2s)));
L= @(theta) (-sum(Ws.*log(F(theta))+(ones(N,1)-W).*log(ones(N,1)-F(theta))));
theta0s=[1;7];
thetahats= fminsearch(L,theta0s);
pXhats=F(thetahats); % Estimated Propensity Score


% first, construct estimate for ATE tau based on estimated propensity score
Ybar1s=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
for i=1:N
    if Ws(i)==1
        DPS= abs( pXhats(i)*W0s-pXhats.*W0s );
        DPS(Ws==1)=2;
        [mi,ind]=min(DPS); 
        Ybar1s(i)=Ys(ind);
    else 
        DPS= abs(pXhats(i)*Ws-pXhats.*Ws);
        DPS(Ws==0)=2;
        [mi,ind]=min(DPS);
        Ybar1s(i)=Ys(ind);
    end
end
tauhats= sum( ( 2*Ws-ones(N,1) ).*(Ys-Ybar1s) )/N; % tauhat based on estimated propensity score 

Tdistr(b)=sqrt(N)*(tauhats-tauhat);
    
end
Crit = invquantile(Tdistr, 1-alpha);

teststat=sqrt(N)*(tauhat-tau);
 if teststat>Crit
     Rej=Rej+1;
 end
 
 
end
RejProb= Rej/R;

fprintf('DGP3 with N= %d, rejection probabity is %1.3f and time is %4.2f mins\n\n',N,  RejProb, time);

% For M=2
% When N=100, R=1000, B=500, the rejection probability from this bootstrap method
% is 0.1270  for DGP3

% When N=200, R=1000, B=500, the rejection probability from this bootstrap method
% is 0.1350 for DGP3

% When N=500, R=1000, B=500, the rejection probability from this bootstrap method
% is 0.090 for DGP3


% When N=500, R=500, B=500, the rejection probability from this bootstrap method
% is 0.078 for DGP3

% When N=1000, R=1000, B=500, the rejection probability from this bootstrap method
% is 0.069 for DGP3

% For M=0
% When N=100, R=1000, B=500, the rejection probability from this bootstrap method
% is   for DGP3

% When N=200, R=1000, B=500, the rejection probability from this bootstrap method
% is  for DGP3

% When N=500, R=1000, B=500, the rejection probability from this bootstrap method
% is  for DGP3


% When N=500, R=500, B=500, the rejection probability from this bootstrap method
% is  for DGP3

% When N=1000, R=1000, B=500, the rejection probability from this bootstrap method
% is  for DGP3

% When N=2000, R=200, B=500, the rejection probability from this bootstrap method
% is 0.065 for DGP3


% Conclusion: This bootstrap method works. 

% Change M to do more experiments. 


