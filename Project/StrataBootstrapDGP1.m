% Strata Bootstrap
% draw X, then draw W using thetahat, then draw Y based on matching. 

%%%% DGP1
%%%% H0: tau<=5 vs H1: tau>5
clear;
rng(246);
N=2000;
alpha=0.05;
B=500;
R=100;
Rej=0;
M=1;

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
theta0=[1;2];
thetahat= fminsearch(L,theta0);
pXhat=F(thetahat); % Estimated Propensity Score

     % compute J_{w}(i)
JNN=zeros(N,1); % compute J_{NN}(i)
 for i=1:N
    if W(i)==0
       Dist= sum(abs( (repmat(X(i,:),N,1)-X).*repmat(W,1,2) ),2);
       Dist(W==0)=inf;
       [~,match]=min(Dist);
       JNN(i)=match;
    else
       Dist= sum(abs( (repmat(X(i,:),N,1)-X).*repmat(W0,1,2) ),2);
       Dist(W0==0)=inf;
       [~,match]=min(Dist);
       JNN(i)=match;
    end
 end
Jw0= W0.*[1:N]'+ W.*JNN; % Jw for w=0
Jw1= W.*[1:N]'+ W0.*JNN; % Jw for w=1
Y0p= Y(Jw0);
Y1p= Y(Jw1);

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


Tdistr=zeros(B,1);
for b=1:B
     S= unidrnd(N,N,1);
     Xb= X(S,:);
     pXb= exp(Xb*thetahat)./( ones(N,1)+ exp(Xb*thetahat) );
     Wb= binornd(1,pXb);
     W0b= ones(N,1)- Wb;
     Yb=Y(S);
     Yb= W0b.*Y0p+ Wb.*Y1p;

% Estimate theta
F= @(theta) (exp(theta(1)*Xb(:,1)+theta(2)*Xb(:,2))./(ones(N,1)+exp(theta(1)*Xb(:,1)+theta(2)*Xb(:,2))));
L= @(theta) (-sum(W.*log(F(theta))+(ones(N,1)-W).*log(ones(N,1)-F(theta))));
theta0=[1;2];
thetahatb= fminsearch(L,theta0);
pXhatb=F(thetahatb); % Estimated Propensity Score

% first, construct estimate for ATE tau based on estimated propensity score
Ybar1b=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if Wb(i)==1
        DPS= abs( pXhatb(i)*W0b-pXhatb.*W0b );
        DPS(Wb==1)=3;
        [mi,ind]=min(DPS); 
        Ybar1b(i)=Yb(ind);
    else 
        DPS= abs(pXhatb(i)*Wb-pXhatb.*Wb);
        DPS(Wb==0)=3;
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


% When N=1000, R=100, B=500, rejection probability using this bootstrap
% method is 0.1000

% When N=500, R=100, B=500, rejection probability using this bootstrap
% method is 0.1600

% when N=1000, R=500, B=500, rejection probability is 0.0780
