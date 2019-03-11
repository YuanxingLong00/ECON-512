%% Empirical Project Adusumilli's JMP 
%% First, Four SImulations with different DGP
% First, Simulation with DGP taken from Abadie and Imbens (2016).
% DGP 1
clear;
rng(115)
N=100;
Rej=0;
alpha=0.05;
Crit= norminv(1-alpha/2);

for k=1:1000
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
K= @(u) ( 0.75*(ones(N,1)-u.^2).*indicator(u) ); % The kernel to be used 
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

% Calculate I
s2= zeros(2,2);
for i=1:N
    Xthe= X(i,:)*thetahat;
    s2=s2+ (f(Xthe)^2 /(pXhat(i)*(1-pXhat(i))))* (X(i,:)'*X(i,:));
end
I= s2/N;
vartau=sigmasq-c*(I\c'); % estimated variance based on estimated propensity score. 

asystat=sqrt(N)*abs(tautilde-5)/sqrt(vartau); % test statistic for two-sided test, alpha=5%
if asystat>Crit
    Rej=Rej+1;
end
end

RejProb=Rej/1000;


%% The simulated rejection probability is 0.02 which is too small.
% Choice of bandwidth matters a lot in this situation.
% need to check!
% The key problem is that kernel estimation does not work in this case, why
% ? Therefore, I go back to series estimation by assuming that E(Y|X)=m(X).


