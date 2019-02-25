%% Empirical Project Adusumilli's JMP 
%% First, Four SImulations with different DGP
% First, Simulation with DGP taken from Abadie and Imbens (2016).
% DGP 1
clear;
rng(135)
N=100;
X1= rand(N,1)-0.5*ones(N,1);
X2= rand(N,1)-0.5*ones(N,1);
X=[X1 X2];
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= 3*X1-3*X2+U0;
Y1= 5+5*X1+X2+U1;
pX= (exp(X1+2*X2))./(ones(N,1)+exp(X1+2*X2));
W= binornd(1,pX);
Y= W.*Y1+(ones(N,1)-W).*Y0;

% Estimate theta
F= @(theta) (exp(theta(1)*X1+theta(2)*X2)./(ones(N,1)+exp(theta(1)*X1+theta(2)*X2)));
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
        DPS= abs( pX(i)*(ones(N,1)-W)-pX.*(ones(N,1)-W) );
        DPS(DPS==0)=2;
        [mi,ind]=min(DPS); 
        Ybar(i)=Y(ind);
    else 
         DPS= abs(pX(i)*W-pX.*W);
         DPS(DPS==0)=2;
        [mi,ind]=min(DPS);
        Ybar(i)=Y(ind);
    end
end
tauhat= sum( (2*W-ones(N,1)).*(Y-Ybar./M) )/N; % tauhat based on known propensity score 


%% first, construct estimate for ATE tau based on estimated propensity score
Ybar1=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if W(i)==1
        DPS= abs( pXhat(i)*(ones(N,1)-W)-pXhat.*(ones(N,1)-W) );
        DPS(DPS==0)=2;
        [mi,ind]=min(DPS); 
        Ybar1(i)=Y(ind);
    else 
        DPS= abs(pXhat(i)*W-pXhat.*W);
        DPS(DPS==0)=2;
        [mi,ind]=min(DPS);
        Ybar1(i)=Y(ind);
    end
end
tautilde= sum( (2*W-ones(N,1)).*(Y-Ybar1) )/N; % tautilde based on estimated propensity score 

%%%% Estimate asymptotic variance of tauhat when propensity score is known but based on estimated prop score 
K= @(u) ( (3*(ones(N,1)-u.^2)/4).*indicator(u) ); % The kernel to be used 
h= 1/sqrt(N); % bandwidth 
mu1= @(p) ( sum( (Y.*W).*K((pXhat-p*ones(N,1))/h) ) / sum( W.*K((pXhat-p*ones(N,1))/h) ) ); % mu(1,p)
mu0= @(p) ( sum( (Y.*(ones(N,1)-W)).*K((pXhat-p*ones(N,1))/h) ) / sum( (ones(N,1)-W).*K((pXhat-p*ones(N,1))/h) ) ); % mu(0,p)
Y2=Y.^2; 
semom1= @(p) ( sum( (Y2.*W).*K((pXhat-p*ones(N,1))/h) ) / sum( W.*K((pXhat-p*ones(N,1))/h) ) ); % E(Y^2|W=1,p(X)=p)
semom0= @(p) ( sum( (Y2.*(ones(N,1)-W)).*K((pXhat-p*ones(N,1))/h) ) / sum( (ones(N,1)-W).*K((pXhat-p*ones(N,1))/h) ) ); % E(Y^2|W=0,p(X)=p)
sigma1= @(p) (semom1(p)-mu1(p)^2); % var(Y|W=1,p(X)=p)
sigma0= @(p) (semom0(p)-mu0(p)^2); % var(Y|W=1,p(X)=p)
s=0;
for i=1:N
    s=s+(mu1(pXhat(i))-mu0(pXhat(i))-tautilde)^2+...
        sigma1(pXhat(i))*(1/pXhat(i)+(1/pXhat(i)-pXhat(i))/(2*M)) + ...
        sigma0(pXhat(i))*( 1/(1-pXhat(i))+(1/(1-pXhat(i))-1+pXhat(i))/(2*M) );
end
sigmasq= s/N;

%%%% estimate the reduction c'Ic in variance when using estimated propensity
%%%% score
f= @(y) ( exp(y)/(1+exp(y))^2);

%%% calculate c
mu1x= @(x) ( sum(  Y.*W.*K((X1-x(1)*ones(N,1))/h).*K( (X2-x(2)*ones(N,1))/h )  )...
           / sum( W.*K((X1-x(1)*ones(N,1))/h).*K((X2-x(2)*ones(N,1))/h) ) );  %% this function generates NaN for some observations!!
mu0x= @(x) ( sum(  Y.*(ones(N,1)-W).*K((X1-x(1)*ones(N,1))/h).*K( (X2-x(2)*ones(N,1))/h)  )...
          /sum( (ones(N,1)-W).*K((X1-x(1)*ones(N,1))/h).*K((X2-x(2)*ones(N,1))/h) ) ); 
mu1xvec= zeros(N,1); % mu(1,X) as a vector
mu0xvec= zeros(N,1); % mu(0,X) as a vector 
for i=1:N
    mu1xvec(i)= mu1x(X(i,:));
    mu0xvec(i)= mu0x(X(i,:));
end

% calculate cov[X,mu(1,X)|F(X'thetatilde)] and cov[X,mu(0,X)|F(X'thetatilde)]
mu1xmat=repmat(mu1xvec,1,2);
mu0xmat=repmat(mu0xvec,1,2);
Xmu1= X.*mu1xmat;
Xmu0= X.*mu0xmat;

% cov[X,mu(1,X)|p] and cov[X,mu(0,X)|p]
EX= @(p) (sum( X.*repmat(K((pXhat-p*ones(N,1))/h),1,2), 1 ) / sum(K((pXhat-p*ones(N,1))/h) ) ); %% E(X|p(X)=p)
EXmu1= @(p) (sum( Xmu1.*repmat(K((pXhat-p*ones(N,1))/h),1,2), 1 ) / sum(K((pXhat-p*ones(N,1))/h) ) ); % there is a problem
EXmu0= @(p) (sum( Xmu0.*repmat(K((pXhat-p*ones(N,1))/h),1,2), 1 ) / sum(K((pXhat-p*ones(N,1))/h) ) ); % there is a problem
Emu1x= @(p) (sum( mu1xmat.*repmat(K((pXhat-p*ones(N,1))/h),1,2), 1 ) / sum(K((pXhat-p*ones(N,1))/h) ) );
Emu0x= @(p) (sum( mu0xmat.*repmat(K((pXhat-p*ones(N,1))/h),1,2), 1 ) / sum(K((pXhat-p*ones(N,1))/h) ) );
CovX1= @(p) ( EXmu1(p)-EX(p).*Emu1x(p) );
CovX0= @(p) ( EXmu0(p)-EX(p).*Emu0x(p) );

% Calculate c
s1=[0,0];
for i=1:N
    Xthe= X(i,:)*thetahat;
    s1=s1+ ( CovX1(pXhat(i))/pXhat(i) + CovX0(pXhat(i))/(1-pXhat(i)) ).*f(Xthe);
end
c= s1/N;


