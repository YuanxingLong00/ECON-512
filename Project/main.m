%% Empirical Project Adusumilli's JMP 
%% First, Four SImulations with different DGP
% First, Simulation with DGP taken from Abadie and Imbens (2016).
% DGP 1
clear;
rng(135)
N=100;
X1= rand(N,1);
X2= rand(N,1);
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= 3*X1-3*X2-U0;
Y1= 5+5+5*X1+X2+U1;
pX= (exp(X1+2*X2))./(1+exp(X1+2*X2));
W= binornd(1,pX);
Y= W.*Y1+(ones(N,1)-W).*Y0;

% Estimate theta
F= @(theta) (exp(theta(1)*X1+theta(2)*X2)./(1+exp(theta(1)*X1+theta(2)*X2)));
L= @(theta) (-sum(W.*log(F(theta))+(ones(N,1)-W).*log(ones(N,1)-F(theta))));
theta0=[1;1];
thetahat= fminsearch(L,theta0);
pXhat=F(thetahat); % Estimated Propensity Score

%%%% Asymptotic method

%% first, construct estimate for ATE tau based on known propensity score
Ybar=zeros(N,1); % i's matcher's average outcome based on known propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if W(i)==1
        [mi,ind]=min(abs(pX(i).*ones(N,1)-pX.*(ones(N,1)-W)) );
        Ybar(i)=Y(ind);
    else 
        [mi,ind]=min( abs(pX(i).*ones(N,1)-pX.*W) );
        Ybar(i)=Y(ind);
    end
end
tauhat= sum( (2*W-ones(N,1)).*(Y-Ybar./M) )/N; % tauhat based on known propensity score 

%% first, construct estimate for ATE tau based on estimated propensity score
Ybar1=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if W(i)==1
        [mi,ind]=min(abs(pXhat(i).*ones(N,1)-pXhat.*(ones(N,1)-W)) );
        Ybar1(i)=Y(ind);
    else 
        [mi,ind]=min( abs(pXhat(i).*ones(N,1)-pXhat.*W) );
        Ybar1(i)=Y(ind);
    end
end
tautilde= sum( (2*W-ones(N,1)).*(Y-Ybar1./M) )/N; % tautilde based on estimated propensity score 

%%%% Estimate asymptotic variance of tauhat when propensity score is known. 
K= @(u) ( (3*(ones(N,1)-u.^2)/4).*indicator(u) ); % The kernel to be used 
h= 1/sqrt(N); % bandwidth 
mu1= @(p) ( sum( (Y.*(1-W)).*K((pX-p*ones(N,1))/h) ) / sum( W.*K((pX-p*ones(N,1))/h) ) ); % mu(1,p)
mu0= @(p) ( sum( (Y.*(1-W)).*K((pX-p*ones(N,1))/h) ) / sum( (1-W).*K((pX-p*ones(N,1))/h) ) ); % mu(0,p)
Y2=Y.^2; 
semom1= @(p) ( sum( (Y2.*(1-W)).*K((pX-p*ones(N,1))/h) ) / sum( W.*K((pX-p*ones(N,1))/h) ) ); % E(Y^2|W=1,p(X)=p)
semom0= @(p) ( sum( (Y2.*(1-W)).*K((pX-p*ones(N,1))/h) ) / sum( (1-W).*K((pX-p*ones(N,1))/h) ) ); % E(Y^2|W=0,p(X)=p)
sigma1= @(p) (semom1(p)-mu1(p)^2); % var(Y|W=1,p(X)=p)
sigma0= @(p) (semom0(p)-mu0(p)^2); % var(Y|W=1,p(X)=p)
s=0;
for i=1:N
    s=s+(mu1(pX(i))-mu0(pX(i))-tauhat)^2+ sigma1(pX(i))*(1/pX(i)+(1/pX(i)-pX(i))/(2*M)) + sigma0(pX(i))*( 1/(1-pX(i))+(1/(1-pX(i))-1+pX(i))/(2*M) );
end
sigmasq= s/N;

