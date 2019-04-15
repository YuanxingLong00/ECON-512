%%%% Bootstrap method for propensity score matching estimator with potential outcomes.
%%%% DGP1
%%%% H0: tau<=5 vs H1: tau>5
clear;
rng(200);
N=1000;
alpha=0.05;
B=1000;
R=2000;
Rej=0;

for r=1:R
% generate data
X1= rand(N,1)-0.5*ones(N,1);
X2= rand(N,1)-0.5*ones(N,1);
X=[X1 X2];
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= 3*X1-3*X2+U0;
Y1= 5+5*X1+X2+U1;
tau=5;
pX= (exp(X1+2*X2))./(ones(N,1)+exp(X1+2*X2));
W= binornd(1,pX);
W0=ones(N,1)-W;
Y= W.*Y1+(ones(N,1)-W).*Y0;
% Estimate theta
F= @(theta) (exp(theta(1)*X1+theta(2)*X2)./(ones(N,1)+exp(theta(1)*X1+theta(2)*X2)));
L= @(theta) (-sum(W.*log(F(theta))+(ones(N,1)-W).*log(ones(N,1)-F(theta))));
theta0=[1;1];
thetahat= fminsearch(L,theta0);
pXhat=F(thetahat); % Estimated Propensity Score

% first, construct estimate for ATE tau based on estimated propensity score
Ybar1=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
M=1; % Use M=1 now. Will use other M later.
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

% Estimate potential outcome. 
Y1= W.*Y+W0.*Ybar1;
Y0= W0.*Y+W.*Ybar1;

Tdistr=zeros(B,1);
for i=1:B
    S= unidrnd(N,N,1);
    Y1b=Y1(S);
    Y0b=Y0(S);
    taub= sum(Y1b-Y0b)/N;
    Tdistr(i)=sqrt(N)*(taub-tauhat);
end
Crit = invquantile(Tdistr, 1-alpha);

teststat=sqrt(N)*(tauhat-tau);
 if teststat>Crit
     Rej=Rej+1;
 end
end
 RejProb= Rej/R



