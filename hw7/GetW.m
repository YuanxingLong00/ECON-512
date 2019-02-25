function [W0,W1,W2]= GetW(V)
L=30;
delta=0.03;

W0=zeros(L,L);
W1=zeros(L,L);
W2=zeros(L,L);
for i=1:L
    for j=1:L
        deltai=1-(1-delta)^i;
        deltaj=1-(1-delta)^j;
        W0(i,j)=V(i,j)*(1-deltai)*(1-deltaj)+V(max(i-1,1),j)*deltai*(1-deltaj)+V(i,max(j-1,1))*(1-deltai)*deltaj+V(max(i-1,1),max(j-1,1))*deltai*deltaj;
        W1(i,j)=V(min(i+1,L),j)*(1-deltai)*(1-deltaj) +V(i,j)*deltai*(1-deltaj)+ V(min(i+1,L),max(j-1,1))*(1-deltai)*deltaj+ V(i,max(j-1,1))*deltai*deltaj;
        W2(i,j)=V(i,min(j+1,L))*(1-deltai)*(1-deltaj) +V(max(i-1,1),min(j+1,L))*deltai*(1-deltaj)+ V(i,j)*(1-deltai)*deltaj+ V(max(i-1,1),j)*deltai*deltaj;
    end
end


