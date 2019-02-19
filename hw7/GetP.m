function np=GetP(W0,W1,W2,p)
global L luo kapa l c v delta beta lambda CRIT D0 D1 D2; 
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


    