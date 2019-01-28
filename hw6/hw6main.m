%% HW6 
%% Question 2
clear
J=21;
I=100;
iter =1000;
[p,P,d]=tauchen(J,1,0.5,0.1,3);

%% Solve by VFI. 
delta=0.95;
s=1:1:I; % state for stock
V= ones(I,J); 
spol=zeros(I,J);
VI= zeros(I,J);
for k=1:iter
   for i=1:I
      for j=1:J
          VV=(ones(I,1).*s(i)-s').*p(j)-0.2*((ones(I,1).*s(i)-s').^1.5) + (V*P(:,j)).*delta;
          VVpartial=VV(1:i,:);
          VI(i,j)=max(VVpartial);
      end
   end
   if max(abs(V-VI))<= 10^(-8)
       V=VI;
       fprintf('k= %d \t Convergence achieved \n', k);
       %to find policy function
          for i=1:I
           for j=1:J
            VV=(ones(I,1).*s(i)-s').*p(j) -0.2*((ones(I,1).*s(i)-s').^1.5) + (V*P(:,j)).*delta;
            VVpartial=VV(1:i,:);
            [VI(i,j),sindex]=max(VVpartial);
            spol(i,j)=s(sindex);
           end
          end
      break
   else
       V=VI;
        fprintf('k= %d \t still not convergent \n', k);
   end
end

% plot the value of the firm depending the initial stocks.
plot(s,V(:,8),'DisplayName','p=0.9')
hold on
plot(s,V(:,11),'DisplayName','p=1')
plot(s,V(:,14),'DisplayName','p=1.1')
hold off
lgd = legend;
lgd.NumColumns = 2;

%% plot the policy function depending on current price
plot(p,spol(25,:),'DisplayName','s=25')
hold on
plot(p,spol(50,:),'DisplayName','s=50')
plot(p,spol(75,:),'DisplayName','s=75')
hold off
lgd = legend;
lgd.NumColumns = 2;       

%% simulate the paths and calculate the expected values and confidence intervals
sim=1000;
% intial state
s0=100;
p0=1;
k0=11; % initial state index for price
T=20; 
simpath=zeros(sim,T);
for s=1:sim
        ppathind=simmar(T,k0,p,P);
        simpath(s,1)=spol(s0,k0);
        for t=2:T
            simpath(s,t)=spol(simpath(s,t-1),ppathind(t));
        end
end
expspath=mean(simpath,1);
lower=quantile(simpath,0.05);
upper = quantile(simpath, 0.95);
Tseq=1:20;
plot(Tseq, expspath,Tseq, lower,Tseq, upper)

            
        


        
