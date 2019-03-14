%%%% This code want to implement the Bootstrap Inference Method.
%%%% DGP1
%%%% H0: tau<=5 vs H1: tau>5
clear;
rng(200);
N=100;
N
alpha=0.05;
M=1;
Rep=100;
Rej=0;
RepM=1;
B=400;
qN=5;

tic

 for R=1:Rep
% Generate Data
Blocksize=zeros(qN,2);
while ~all(Blocksize(:)~=0) % Eliminate unbalanced sample 
    
X1= rand(N,1)-0.5*ones(N,1);
X2= rand(N,1)-0.5*ones(N,1);
X=[X1 X2];
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= 3*X1-3*X2+U0;
Y1= 5*ones(N,1)+5*X1+X2+U1;
tau=5;
pX= (exp(X1+2*X2))./(ones(N,1)+exp(X1+2*X2));
W= binornd(1,pX);
W0= ones(N,1)-W;
Y= W.*Y1+W0.*Y0;


% Estimate theta
F= @(theta) (exp(theta(1)*X1+theta(2)*X2)./(ones(N,1)+exp(theta(1)*X1+theta(2)*X2)));
L= @(theta) (-sum(W.*log(F(theta))+(ones(N,1)-W).*log(ones(N,1)-F(theta))));
% fF= @(theta)((ones(N,1)-(exp(theta(1)*X1+theta(2)*X2)))./(ones(N,1)+exp(theta(1)*X1+theta(2)*X2)) );
% matrix= @(theta)( W.*fF(theta)-W0.*fF(theta).*exp(theta(1)*X1+theta(2)*X2) );
% FOC = @(theta)(- sum( (X.* repmat(matrix(theta),1,2)), 1)  );
theta0=[0.9;1.9];
% options = optimset('Display','off');
thetahat= fminsearch(L,theta0);
pXhat=F(thetahat); % Estimated Propensity Score

% This part is to eliminate unbalanced original sample. 
qN=5;
len=floor(N/qN);
[fvalue,xvalue]=ecdf(pXhat);
xval=xvalue(2:(N+1));
PSranking=fvalue(2:(N+1));
Index= zeros(N,1);
iblock=zeros(N,1);
for i=1:N
    Index(i)= find(pXhat==xval(i));
end
Wr= W(Index);
block= zeros((qN+1),1);
for i=1:qN
    block(i+1)= find(PSranking<=(i/qN),1,'last');
end


for i=1:N
    ipsrank= find(Index==i);
    block(1)=0;
   if ipsrank <N
        iblock(i)= find(block>ipsrank,1,'first')-1;
    else
        iblock(i)=qN;
   end
    
    if W(i)==0
        Blocksize(iblock(i),1)=Blocksize(iblock(i),1)+1;
    else 
        Blocksize(iblock(i),2)=Blocksize(iblock(i),2)+1;
    end
end
if ~all(Blocksize(:)~=0)
    Blocksize=zeros(qN,2);
    continue
end
  
  break
end

%% first, construct estimate for ATE tau based on thetahat
Ybar=zeros(N,1); % i's matcher's average outcome based on thetahat
for i=1:N
    if W(i)==1
        DPS= abs( pXhat(i)*W0-pXhat.*W0 );
        DPS(W0==0)=3;
        [mi,ind]=min(DPS); 
        Ybar(i)=Y(ind);
    else 
         DPS= abs(pXhat(i)*W-pXhat.*W);
         DPS(W==0)=3;
        [mi,ind]=min(DPS);
        Ybar(i)=Y(ind);
    end
end
tauhat= sum( (2*W-ones(N,1)).*(Y-Ybar) )/N; % tauhat based on thetahat
% This will be used in inference. 

% Implement Bootstrap Method 

%-----------Step 0: Obtain M(i) and Jw(i)----------------------------%

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
Jw=zeros(N,2);
Jw(:,1)= W0.*[1:N]'+ W.*JNN; % Jw for w=0
Jw(:,2)= W.*[1:N]'+ W0.*JNN; % Jw for w=1


% Obtain Mi which is the multinominal distribution realizations. It is
% important to obatain it in this application.qN=N^(1/3). qN cannot be too
% large.

        
Swli= zeros(N, 3*len); % Each i's match from the other group
Mi= zeros(N, 3*len);   % Standarized multinominal probabilities
Matsiz=zeros(N,1); % Number of i's match 
block(1)=1;
for i=1:N
    ipsrank= find(Index==i);
    if iblock(i) < qN % to be adjusted for other codes. 
         imatchl= block(iblock(i));
         imatchu= block((iblock(i)+1))-1;
    else 
        imatchl= block(iblock(i));
        imatchu= block((iblock(i)+1));
    end
    
    blo_ps= Wr( imatchl: imatchu );
    blo_ind= Index( imatchl: imatchu );
    blo_ps0= 1-blo_ps;
    if W(i)==0
       siz= sum(blo_ps);
       temp= (blo_ps.* blo_ind)';
       Swli(i,1:siz)= temp(~(temp==0));
       pmn= ones(1,siz)/siz;
       temp=mnrnd(siz,pmn);
       Mi(i,1:siz)= temp./siz;
       Matsiz(i)=siz;
    else 
       siz= sum(blo_ps0);
       temp= (blo_ps0.* blo_ind)';
       Swli(i,1:siz)= temp(~(temp==0));
       pmn=ones(1,siz)/siz;
       temp=mnrnd(siz,pmn);
       Mi(i,1:siz)= temp./siz;
       Matsiz(i)=siz;
    end
end


        



Tdistr=zeros(B,1);
for sim=1:B
%-----------Step 1: Sample covariates Xs  ----------------------------%
N1s= (M+1);
N0s= (M+1);

while N1s<= (M+1) || N0s<=(M+1) % Eliminate unbalanced bootstrap sample 
S= unidrnd(N,N,1);
Xs =zeros(N,2);
for i=1:N
    Xs(i,:) = X(S(i),:);
end

%-----------Step 2: Draw new treatment values Ws----------------------%
pXhats= exp(thetahat(1)*Xs(:,1)+thetahat(2)*Xs(:,2))./(ones(N,1)+exp(thetahat(1)*Xs(:,1)+thetahat(2)*Xs(:,2)));
Ws= binornd(1,pXhats);
W0s= ones(N,1)-Ws;
%-----------Step 3: Discard unbalanced samples------------------------%
 N1s= sum(Ws);
 N0s= N- N1s;
 if  N1s<= (M+1) || N0s<=(M+1)
     continue
 end
 break
end


% Estimate theta based bootstrap sample
F= @(theta) (exp(theta(1)*Xs(:,1)+theta(2)*Xs(:,2))./(ones(N,1)+exp(theta(1)*Xs(:,1)+theta(2)*Xs(:,2))));
L= @(theta) (-sum(Ws.*log(F(theta))+W0s.*log(ones(N,1)-F(theta))));
theta0=[0.9;1.9];
thetahats= fminsearch(L,theta0);
% pXhats=F(thetahats); % Estimated Propensity Score

%-----------Step 4: Calculate the multiple objects-------------------%
% compute K_M(i,thetahats)
pXhattemp= exp(thetahats(1)*X(:,1)+thetahats(2)*X(:,2))./(ones(N,1)+exp(thetahats(1)*X(:,1)+thetahats(2)*X(:,2)));
Matchset=zeros(N,1); % Match from opposite group for any i=1,...N
for i=1:N
    if W(i)==1
        DPS= abs( pXhattemp(i)*W0-pXhattemp.*W0 );
        DPS(W0==0)=2;
        [mi,ind]=min(DPS); 
        Matchset(i)=ind;
    else 
        DPS= abs(pXhattemp(i)*W-pXhattemp.*W);
        DPS(W==0)=2;
        [mi,ind]=min(DPS);
        Matchset(i)=ind;
    end
end
KM=zeros(N,1); % It is K_M(i,thetahats)
for i=1:N
    KM(i)= sum(Matchset(:) == i);
end


% compute KMtilde where w \neq W_i  (Comment: this part is currently
% wrong.)
KMtilde=zeros(N,1);
for i=1:N
    imatchvec= Swli(i,:)';
    imatch= imatchvec(~(imatchvec==0));
    imatchsiz= size(imatch);
    imatchprob= Mi(i,1:imatchsiz);
    KMtilde(i)= KM(imatch)'* imatchprob';
end

KMtil= zeros(N,2); % it contains all the matching functions to be used. 
KMtil(:,1)= W0.*KM+ W.*KMtilde; % KMtil for w=0
KMtil(:,2)= W.*KM+ W0.*KMtilde; % KMtil for w=1






% compute tauhats
%% first, construct estimate for ATE tau based on thetahats
Ybar=zeros(N,1); % i's matcher's average outcome based on thetahats
for i=1:N
    if W(i)==1
        DPS= abs( pXhattemp(i)*W0-pXhattemp.*W0 );
        DPS(W0==0)=3;
        [mi,ind]=min(DPS); 
        Ybar(i)=Y(ind);
    else 
         DPS= abs(pXhattemp(i)*W-pXhattemp.*W);
         DPS(W==0)=3;
        [mi,ind]=min(DPS);
        Ybar(i)=Y(ind);
    end
end
tauhats= sum( (2*W-ones(N,1)).*(Y-Ybar./M) )/N; % tauhat based on thetahats


% compute e_{2i}(w,thetahats) and compute e_1(i)
% first compute e_{2i}(W_i, thetahats) and e_{2i}(thetahats)
e2Wi=zeros(N,1);
e1=zeros(N,1);

%% series estimation linear for muwip0 and muwip1
serieslinear;

%% series estimation 3 polinominals for muwip0 and muwip1
% series3polinominals;
%% series estimation 4 polinominals for muwip0 and muwip1
% series3polinominals;

%% Kernel estimation for muwip0 and muwip1
% serieskernel; 


e2=zeros(N,2);

e2(:,1)=e2Wi(Jw(:,1));
e2(:,2)=e2Wi(Jw(:,2));

% compute nu
nu= (1+KMtil/M).*e2;

% compute bootstrap realized error

nu0= nu(:,1);
nu1= nu(:,2); 
bterror= e1(S)+Ws.*nu1(S)- W0s.* nu0(S);


% compute the center of bootstrap realized error 
CounterpX= ones(N,1)- pXhattemp;
center= e1+ pXhattemp.*nu1- CounterpX.* nu0;

%-----------------Step 5: compute the test statistic--------------%
t = sum(bterror- center)/sqrt(N);
Tdistr(sim)=t;
end
Crit = invquantile(Tdistr, 1-alpha);

teststat=sqrt(N)*(tauhat-tau);
 if teststat>Crit
     Rej=Rej+1;
 end
 
 if R==250
     R
 else 
     if R==500
         R
     else if R==750
             R
         end
     end
 end
     
 end
 RejProb= Rej/Rep
 toc 
 time=toc/60

% The key idea of this bootstrap method is to find the critical value for
% tauhat and use the fact that statistic T has the same asymptotic
% distribution as sqrt(N)(tauhat-tau). In other words, the bootstrap
% procedure is consistent. 
% Rep=100 might be too small. Need to check for a larger Simulation number.

% How to deal with unbalanced sample (W,X) and (W^*,X^*)? Drop it?

% Rejction Probability when N=100 is 0.087 when using kernel qN=5
% time = 30.3472 mins
% Rejction Probability when N=200 is 0.076 when using kernel 
% time = 3.4637e+03 sec
% Rejction Probability when N=500 is 
% Rejction Probability when N=1000 is 


% Rejction Probability when N=100 is 0.5310 when using series linear qN=5
% time = 20.88 mins
% Rejction Probability when N=200 is  when using kernel linear qN=5
% time = 

% Rejction Probability when N=500 is 
% Rejction Probability when N=1000 is 



