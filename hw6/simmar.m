function path = simmar(T,k0,p,P)
% N    the length of simulated markov process
% k0   initial state (index)
% p    the list of states
% P    transition matrix
shock= rand(T,1);
% initial values of the indicator
zz=zeros(T,1);
% suppose R start from state 2 (so the mean value of R)
zz(1)=k0;
Pcumu=cumsum(P',2);
for t=2:T
    if shock(t)<Pcumu(zz(t-1),1)
        zz(t)=1;
    else
        k=find(Pcumu(zz(t-1),:)<shock(t),1,'last');
        zz(t)=k+1;
    end
end
path=p(zz);




