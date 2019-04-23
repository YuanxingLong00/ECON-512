%%%% Propensity score matching 
%%%% DGP0 
%%%% H0: mu<=5 VS H1: mu>5
clear;
rng(200);
N=100;
alpha=0.05;
M=1;
qN=5;
Rep=1000;
B=400;
Rej=0;



 for R=1:Rep
% Generate Data
Blocksize=zeros(qN,2);
while ~all(Blocksize(:)~=0) % Eliminate unbalanced sample 

X= rand(N,1)-0.5*ones(N,1);
U0= normrnd(0,1,N,1);
U1= normrnd(0,1,N,1);
Y0= -3*X+U0;
Y1= 5+7*X+U1;
theta=2;
tau=5;
pX= (exp(X*theta))./(ones(N,1)+exp(X*theta));
W= binornd(1,pX);
W0= ones(N,1)-W;
Y= W.*Y1+W0.*Y0;



% Estimate theta
F= @(theta) ( exp(X*theta)./(ones(N,1)+exp(X*theta)) );
L= @(theta) (-sum(W.*log(F(theta))+W0.*log(ones(N,1)-F(theta))));
% fF= @(theta)((ones(N,1)-(exp(theta(1)*X1+theta(2)*X2)))./(ones(N,1)+exp(theta(1)*X1+theta(2)*X2)) );
% matrix= @(theta)( W.*fF(theta)-W0.*fF(theta).*exp(theta(1)*X1+theta(2)*X2) );
% FOC = @(theta)(- sum( (X.* repmat(matrix(theta),1,2)), 1)  );
theta0=4;
% options = optimset('Display','off');
thetahat= fminsearch(L,theta0);
pXhat=F(thetahat); % Estimated Propensity Score


% This part is to eliminate unbalanced original sample. 
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
       Dist= sum(abs( (repmat(X(i),N,1)-X).*W ),2);
       Dist(W==0)=inf;
       [~,match]=min(Dist);
       JNN(i)=match;
    else
       Dist= sum(abs( (repmat(X(i),N,1)-X).*W0 ),2);
       Dist(W0==0)=inf;
       [~,match]=min(Dist);
       JNN(i)=match;
    end
 end
Jw=zeros(N,2);
Jw(:,1)= W0.*[1:N]'+ W.*JNN; % Jw for w=0
Jw(:,2)= W.*[1:N]'+ W0.*JNN; % Jw for w=1


       
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
Xs= X(S);
%-----------Step 2: Draw new treatment values Ws----------------------%
pXhats= exp(Xs*thetahat)./(ones(N,1)+exp(Xs*thetahat));
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
F= @(theta) (exp(Xs*theta)./(ones(N,1)+exp(Xs*theta)));
L= @(theta) (-sum(Ws.*log(F(theta))+W0s.*log(ones(N,1)-F(theta))));
theta0=4;
thetahats= fminsearch(L,theta0);
% pXhats=F(thetahats); % Estimated Propensity Score

%-----------Step 4: Calculate the multiple objects-------------------%
% compute K_M(i,thetahats)
pXhattemp= exp(X*thetahats)./(ones(N,1)+exp(X*thetahats));
Matchset=zeros(N,1); % Match from opposite group for any i=1,...N
for i=1:N
    if W(i)==1
        DPS= abs( pXhattemp(i)*W0-pXhattemp.*W0 );
        DPS(W==1)=2;
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
% serieslinear;

%% series estimation 3 polinominals for muwip0 and muwip1
% series3polinominals;
%% series estimation 4 polinominals for muwip0 and muwip1
% series3polinominals;

%% Kernel estimation for muwip0 and muwip1
serieskernel; 


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
center= mean(e1+ pXhattemp.*nu1- CounterpX.* nu0);

%-----------------Step 5: compute the test statistic--------------%
t = sum(bterror- center)/sqrt(N);
Tdistr(sim)=t;
end
Crit = invquantile(Tdistr, 1-alpha);

teststat=sqrt(N)*(tauhat-tau);
 if teststat>Crit
     Rej=Rej+1;
 end
     
 end
 RejProb= Rej/Rep;
 toc 
 time=toc/60;
fprintf('DGP0 with N= %d and qN=%d, rejection probabity is %1.3f and time is %4.2f mins\n\n',N, qN, RejProb, time);

      









