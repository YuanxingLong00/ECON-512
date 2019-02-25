function np=GetP(W0,W1,W2,p)

L=30;
luo=0.85;
kapa=10;
l=15;
beta=1/1.05;
v=10;

% this tolerance is pretty tight, takes long time to converge. for this you
% might need to optimize speed, not using global variables. 

c=zeros(L,1);
c(1:l)=kapa.*([1:l]'.^(log(luo)/ log(2)) );
c(l+1:L)= c(l).*ones(L-l,1);
D0= @(p1,p2)( 1 / (1 + exp(v-p1) +exp(v-p2)) );
D1= @(p1,p2)( exp(v-p1) / (1 + exp(v-p1) +exp(v-p2)) );
D2= @(p1,p2)( exp(v-p2) / (1 + exp(v-p1) +exp(v-p2)) );


np=zeros(L,L);
for w1=1:L
    for w2=1:L
        % f=@(p1)(  -D1(p1,p(w2,w1))*(p1-c(w1))- beta*( D0(p1,p(w2,w1))*W0(w1,w2)+D1(p1,p(w2,w1))*W1(w1,w2)+ D2(p1,p(w2,w1))*W2(w1,w2) ) ); 
        f=@(p1)( 1-(1-D1(p1, p(w2,w1)))*(p1-c(w1))-beta*W1(w1,w2)+beta*( D0(p1,p(w2,w1))*W0(w1,w2)+ D1(p1,p(w2,w1))*W1(w1,w2)+D2(p1,p(w2,w1))*W2(w1,w2)) ) ;
        p0=10;
        options = optimset('Display','off');
        np(w1,w2)=fsolve(f,p0,options);
    end
end
end


    