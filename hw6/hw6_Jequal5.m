%% HW6 
%% Question 2
clear
J=5;
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
          VV=(ones(I,1).*s(i)-s').*p(j) + (V*P(:,j)).*delta;
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
            VV=(ones(I,1).*s(i)-s').*p(j) + (V*P(:,j)).*delta;
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
plot(s,V(:,3),'DisplayName','p=0.82')
hold on
plot(s,V(:,4),'DisplayName','p=1')
plot(s,V(:,5),'DisplayName','p=1.17')
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
