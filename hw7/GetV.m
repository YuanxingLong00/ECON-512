function nV=GetV(p,np,W0,W1,W2)
global L luo kapa l c v delta beta lambda CRIT D0 D1 D2; 
nV=zeros(L,L);
for w1=1:L
    for w2=1:L
        nV(w1,w2)= D1(np(w1,w2),p(w2,w1))*(np(w1,w2)-c(w1))+beta*( D0(np(w1,w2),p(w2,w1))*W0(w1,w2)+D1(np(w1,w2),p(w2,w1))*W1(w1,w2)+ D2(np(w1,w2),p(w2,w1))*W2(w1,w2) );
    end
end
