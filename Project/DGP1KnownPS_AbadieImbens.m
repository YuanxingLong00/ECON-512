%% This exercise wants to replicate the results of Abadie and Imbens (2016) when the
%% propensity score is known and uses the asymptotic theorem. 
%% It's an inference exercise!
clear;
rng(135)
N=100;
rej=0;
alpha=0.05;
Crit= norminv(1-alpha/2);
for j=1:1000
X1= rand(N,1)-0.5*ones(N,1);
X2= rand(N,1)-0.5*ones(N,1);
X=[X1 X2];
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= 3*X1-3*X2+U0;
Y1= 5+5*X1+X2+U1;
pX= (exp(X1+2*X2))./(ones(N,1)+exp(X1+2*X2));
W= binornd(1,pX);
W0=ones(N,1)-W;
Y= W.*Y1+(ones(N,1)-W).*Y0;

%% first, construct estimate for ATE tau based on known propensity score
Ybar=zeros(N,1); % i's matcher's average outcome based on known propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if W(i)==1
        DPS= abs( pX(i)*W0-pX.*W0 );
        DPS(W0==0)=2;
        [mi,ind]=min(DPS); 
        Ybar(i)=Y(ind);
    else 
        DPS= abs(pX(i)*W-pX.*W);
        DPS(W==0)=2;
        [mi,ind]=min(DPS);
        Ybar(i)=Y(ind);
    end
end
tauhat= sum( (2*W-ones(N,1)).*(Y-Ybar) )/N; % tauhat based on known propensity score 


%%%% Estimate asymptotic variance of tauhat when propensity score is known but based on estimated prop score 
K= @(u) ( (3*(ones(N,1)-u.^2)/4).*indicator(u) ); % The kernel to be used 
h= 1/sqrt(N); % bandwidth 
mu1= @(p) ( sum( (Y.*W).*K((pX-p*ones(N,1))/h) ) / sum( W.*K((pX-p*ones(N,1))/h) ) ); % mu(1,p)
mu0= @(p) ( sum( (Y.*(1-W)).*K((pX-p*ones(N,1))/h) ) / sum( (1-W).*K((pX-p*ones(N,1))/h) ) ); % mu(0,p)
Y2=Y.^2; 
semom1= @(p) ( sum( (Y2.*W).*K((pX-p*ones(N,1))/h) ) / sum( W.*K((pX-p*ones(N,1))/h) ) ); % E(Y^2|W=1,p(X)=p)
semom0= @(p) ( sum( (Y2.*(1-W)).*K((pX-p*ones(N,1))/h) ) / sum( (1-W).*K((pX-p*ones(N,1))/h) ) ); % E(Y^2|W=0,p(X)=p)
sigma1= @(p) (semom1(p)-mu1(p)^2); % var(Y|W=1,p(X)=p)
sigma0= @(p) (semom0(p)-mu0(p)^2); % var(Y|W=1,p(X)=p)
s=0;
for i=1:N
    s=s+(mu1(pX(i))-mu0(pX(i))-tauhat)^2+ sigma1(pX(i))*(1/pX(i)+(1/pX(i)-pX(i))/(2*M)) + sigma0(pX(i))*( 1/(1-pX(i))+(1/(1-pX(i))-1+pX(i))/(2*M) );
end
sigmasq= s/N;

Z=sqrt(N)*abs(tauhat-5)/sqrt(sigmasq);
if Z> Crit
    rej=rej+1;
end
end
RejProb=rej/1000;

%% The rejection probability is 4.60%
