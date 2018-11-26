%% HW 4 Calculating pi
% all is good
%% Question 1 using quasi-Monte Carlo
clear
rng(135)
throw=100000;
t=0;
for i=1:throw
  X=rand(2,1);
  if (X(1)^2+X(2)^2) <= 1
      t=t+1;
  end
end
pi=4*t/throw

%% Question 2 using Newton Cotes 
N=1000;
h=1/N;
t=0;
% 
for i=0.5:1:(N-0.5)
  for j=0.5:1:(N-0.5)
      if h*sqrt(i^2+j^2) <=1
          t=t+1;
      end
  end
end
pi= 4*t/N^2
     
    
    
    
    

%% Question 3 
clear
rng(135)
throw=100000;
t=0;
for i=1:throw
  x=rand(1,1);
  t=sqrt(1-x^2)+t;
end
pi= 4*t/throw

%% Question 4
clear
N=100000;
f=@(x) (sqrt(1-x^2));
a=0;
b=1;
pi=4*int_NC(f,a,b,N)

%% Question 5
% Dart-Throwing

throw=1000;
se=0;
for k=1:200
t=0;
for i=1:throw
  X=rand(2,1);
  if (X(1)^2+X(2)^2) <= 1
      t=t+1;
  end
end
pie=4*t/throw;
se= (pie-pi)^2+se;
end
mse=se/200

% Quasi-MC

throw=1000;

rng(135)
se=0;
for k=1:200
t=0;
for i=1:throw
  x=rand(1,1);
  t=sqrt(1-x^2)+t;
end
pi= 4*t/throw
se= (pie-pi)^2+se;
end
mse=se/200

% NC midpoint
clear
N=100000;
se=0;
for i=1:200
f=@(x) (sqrt(1-x^2));
a=0;
b=1;
pie=4*int_NC(f,a,b,N);
se=(pie-pi)^2 +se;
end
mse = se/200

