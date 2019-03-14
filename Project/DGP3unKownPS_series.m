% DGP 3 one-sided test
clear;
rng(115)
N=1000;
N
Rej=0;
alpha=0.05;
Crit= norminv(1-alpha);
tic
for k=1:1000
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
tauhat= sum( (2*W-ones(N,1)).*(Y-Ybar./M) )/N; % tauhat based on known propensity score 


%% first, construct estimate for ATE tau based on estimated propensity score
Ybar1=zeros(N,1); % i's matcher's average outcome based on estimated propensity score
M=1; % Use M=1 now. Will use other M later.
for i=1:N
    if W(i)==1
        DPS= abs( pXhat(i)*W0-pXhat.*W0 );
        DPS(W0==0)=2;
        [mi,ind]=min(DPS); 
        Ybar1(i)=Y(ind);
    else 
        DPS= abs(pXhat(i)*W-pXhat.*W);
        DPS(W==0)=2;
        [mi,ind]=min(DPS);
        Ybar1(i)=Y(ind);
    end
end
tautilde= sum( (2*W-ones(N,1)).*(Y-Ybar1) )/N; % tautilde based on estimated propensity score 

%%%% Estimate asymptotic variance of tauhat when propensity score is known but based on estimated prop score 
K= @(u) ( 0.75*(ones(N,1)-u.^2).*indicator(u) ); % The kernel to be used 
h= 1/sqrt(N); % bandwidth 

mu1= @(p) ( sum( (Y.*W).*K((pXhat-p*ones(N,1))/h) ) ./ sum( W.*K((pXhat-p*ones(N,1))/h) ) ); % mu(1,p)
mu0= @(p) ( sum( (Y.*(ones(N,1)-W)).*K((pXhat-p*ones(N,1))/h) ) ./ sum( (ones(N,1)-W).*K((pXhat-p*ones(N,1))/h) ) ); % mu(0,p)
Y2=Y.^2; 
semom1= @(p) ( sum( (Y2.*W).*K((pXhat-p*ones(N,1))/h) ) ./ sum( W.*K((pXhat-p*ones(N,1))/h) ) ); % E(Y^2|W=1,p(X)=p)
semom0= @(p) ( sum( (Y2.*(ones(N,1)-W)).*K((pXhat-p*ones(N,1))/h) ) ./ sum( (ones(N,1)-W).*K((pXhat-p*ones(N,1))/h) ) ); % E(Y^2|W=0,p(X)=p)
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

% Compute c
% Using linear series estimation to compute mu(1,X) and mu(0,X). The model
% is correctly specified. 
treated=find(W);
control=find(W0);
sizet=size(treated);
sizec=size(control);
Xt=[ones(sizet(1),1) X(treated,:)];
Xc=[ones(sizec(1),1) X(control,:)];
Xexpand= [ones(N,1) X];
XXt=Xt'*Xt;
XXc=Xc'*Xc;
beta1=XXt\Xt'*Y(treated);
beta2=XXc\Xc'*Y(control);
mu1x= Xexpand*beta1; % mu(1,X)
mu0x= Xexpand*beta2; % mu(0,X)

Xmu1xp= @(p)( sum(X.* repmat( (mu1x.*K((pXhat-p*ones(N,1))/h)) ,1,2), 1) ./ sum( K((pXhat-p*ones(N,1))/h) ) );
Xmu0xp= @(p)( sum(X.*repmat( (mu0x.*K((pXhat-p*ones(N,1))/h)) ,1,2), 1) ./ sum( K((pXhat-p*ones(N,1))/h) ) );
Xp= @(p)( sum(X.*repmat( K((pXhat-p*ones(N,1))/h) ,1,2), 1) ./ sum( K((pXhat-p*ones(N,1))/h) ) );
mu1xp= @(p)( sum(mu1x.*K((pXhat-p*ones(N,1))/h) ) / sum( K((pXhat-p*ones(N,1))/h) ) );
mu0xp= @(p)( sum(mu0x.*K((pXhat-p*ones(N,1))/h) ) / sum( K((pXhat-p*ones(N,1))/h) ) );

cov1p= @(p)( Xmu1xp(p)- Xp(p).*repmat(mu1xp(p),1,2) );
cov0p= @(p)( Xmu0xp(p)- Xp(p).*repmat(mu0xp(p),1,2) );
   
% Calculate I and c
s1=zeros(1,2);
s2= zeros(2,2);
for i=1:N
    Xthe= X(i,:)*thetahat;
    s2=s2+ (f(Xthe)^2 /(pXhat(i)*(1-pXhat(i))))* (X(i,:)'*X(i,:));
    s1= s1+ f(Xthe) * ( cov1p(pXhat(i))/pXhat(i)+  cov0p(pXhat(i))/(1-pXhat(i))    );
end
I= s2/N;
cT=s1/N;
vartau=sigmasq-cT*(I\cT'); % estimated variance based on estimated propensity score. 

 asystat=sqrt(N)*(tautilde-tau)/sqrt(vartau); % test statistic for two-sided test, alpha=5%
if asystat>Crit
   Rej=Rej+1;
end
end
toc 
time =toc/60
RejProb=Rej/1000



% The key point is that we need to estimate the conditonal mean by using
% series estimation method instead of using kernels. 

% Rejection Probabilities: 
% If use N=100,  0.079,  time =  0.3227 min
% If use N=200,  0.087,  time =  0.6447 min
% If use N=500,  0.084,  time =  2.1597 mins
% If use N=1000, 0.065, time =   7.3425 mins

