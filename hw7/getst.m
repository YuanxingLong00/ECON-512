function w=getst(w0,N,p)
% get state after N periods
global L luo kapa l c v delta beta lambda CRIT D0 D1 D2; 
w1=w0(1);
w2=w0(2);
for i=1:N
p1=p(w1,w2);
p2=p(w2,w1);
D1v=D1(p1,p2);
D2v=D2(p1,p2);
t=rand(1);
if t<D1v
    q1=1;
    q2=0;
else if D1v<=t && t<(D1v+D2v)
        q1=0;
        q2=1;
    else q1=0;
        q2=0;
    end
end
delta1=1-(1-delta)^w1;
delta2=1-(1-delta)^w2;
t1=rand(1);
t2=rand(1);
if t1<delta1
    w1=max(1,w1+q1-1);
else w1=min(w1+q1,L);
end
if t2<delta2
    w2=max(1,w2+q2-1);
else w2=min(L,w2+q2);
end


end
w=[w1 w2];
        
