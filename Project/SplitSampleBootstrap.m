%%%% Bootstrap method for propensity score matching estimator with potential outcomes.
%%%% DGP1
%%%% H0: tau<=5 vs H1: tau>5
clear;
rng(246);
N=500;
alpha=0.05;
B=500;
R=500;
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
Y1= W.*Y;
Y1= Y1(W==1);
Y0= W0.*Y;
Y0= Y0(W==0);
Xs1= X(W==1,:);
Xs0= X(W==0,:);
s1= size(Y1,1);
s0= size(Y0,1);

Tdistr=zeros(B,1);
for b=1:B
    S1= unidrnd(s1,s1,1);
    S0= unidrnd(s0,s0,1);
    Yb=[Y1(S1)' Y0(S0)']';
    Xs1b= Xs1(S1,:)';
    Xs2b= Xs0(S0,:)';
    Xb=[ Xs1b Xs2b ]';
    Xb1=Xb(:,1);
    Xb2=Xb(:,2);
    Wb=[ones(1,s1),zeros(1,s0)]';
    W0b= ones(N,1)-Wb;
    
    % Estimate theta
F= @(theta) (exp(theta(1)*Xb1+theta(2)*Xb2)./(ones(N,1)+exp(theta(1)*Xb1+theta(2)*Xb2)));
L= @(theta) (-sum(Wb.*log(F(theta))+W0b.*log(ones(N,1)-F(theta))));
theta0b=[1;1];
thetahatb= fminsearch(L,theta0b);
pXhatb=F(thetahatb); % Estimated Propensity Score

% first, construct estimate for ATE tau based on estimated propensity score
Ybar1b=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if Wb(i)==1
        DPS= abs( pXhatb(i)*W0b-pXhatb.*W0b );
        DPS(Wb==1)=2;
        [mi,ind]=min(DPS); 
        Ybar1b(i)=Yb(ind);
    else 
        DPS= abs(pXhatb(i)*Wb-pXhatb.*Wb);
        DPS(Wb==0)=2;
        [mi,ind]=min(DPS);
        Ybar1b(i)=Yb(ind);
    end
end
tauhatb= sum( ( 2*Wb-ones(N,1) ).*(Yb-Ybar1b) )/N; % tauhatb based on estimated propensity score 

Tdistr(b)=sqrt(N)*(tauhatb-tauhat);
end
Crit = invquantile(Tdistr, 1-alpha);

teststat=sqrt(N)*(tauhat-tau);
 if teststat>Crit
     Rej=Rej+1;
 end
 r
 
end
RejProb= Rej/R


% when N=500, R=500, rejection probability from asymptotic method is 0.048
% and that from this bootstrap method is 0.104 when B=500.
% Split Sample Boostrap has smaller variance in tauhatb. 

% This bootstrap is called naive bootstrap in the Abadie and Imbens (2008).
% This bootstrap does not work for sure. 





