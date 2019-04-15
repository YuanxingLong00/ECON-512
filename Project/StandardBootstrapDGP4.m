%%%% Standard Bootstrap method for propensity score matching estimator.
%%%% DGP4
%%%% H0: tau<=5 vs H1: tau>5
% clear
rng(246);
% N=100
alpha=0.05;
B=500;
R=10;
Rej=0;
M=0;

for r=1:R
% generate data
N1=M;
N0=M;
while N1<= M || N0<=M
X=normrnd(0,1,N,4);
X1=X(:,1);
X2=X(:,2);
X3=X(:,3);
X4=X(:,4);
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= U0;
Y1= 210+27.4*X(:,1)+13.7*X(:,2)+13.7*X(:,3)+13.7*X(:,4) +U1;
tau=210;
theta=[-1;0.5;-0.25;-0.1];
pX= (exp(X*theta))./(ones(N,1)+exp(X*theta) );
W= binornd(1,pX);
W0= ones(N,1)-W;
N1=sum(W);
N0=N-N1;
if  N1<= M || N0<=M
     continue
 end
 break
end

Y= W.*Y1+W0.*Y0;

% Estimate theta
F= @(theta) (exp(theta(1)*X1+theta(2)*X2+theta(3)*X3+theta(4)*X4)./...
    (ones(N,1)+exp(theta(1)*X1+theta(2)*X2+theta(3)*X3+theta(4)*X4)));
L= @(theta) (-sum(W.*log(F(theta))+W0.*log(ones(N,1)-F(theta))));
theta0=[-1;0.5;-0.25;-0.1];
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
X3s= Xs(:,3);
X4s= Xs(:,4);

% Estimate theta
F= @(theta) (exp(theta(1)*X1s+theta(2)*X2s+theta(3)*X3s+theta(4)*X4s)./...
    (ones(N,1)+exp(theta(1)*X1s+theta(2)*X2s+theta(3)*X3s+theta(4)*X4s)));
L= @(theta) (-sum(W.*log(F(theta))+W0.*log(ones(N,1)-F(theta))));
theta0s=[-1;0.5;-0.25;-0.1];
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

fprintf('DGP4 with N= %d, rejection probabity is %1.3f and time is %4.2f mins\n\n',N, RejProb, time);
