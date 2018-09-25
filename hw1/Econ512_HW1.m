%%%% Econ 512 Homework 1 
%% Question 1
X=[1 1.5 3 4 5 7 9 10];
Y1=-2+0.5*X;
Y2=-2+0.5*X.^2;
figure
plot(X,Y1,'*',X,Y2,'*')

%% Question 2
step=(20-(-10))/199;
X_2=-10:step:20;
s=sum(X_2)


%% Question 3
A=[2, 4, 6; 
    1, 7, 5; 
    3, 12, 4];
b=[-2; 3; 10];
C=A'*b
D=(A'*A)\b
E=sum(C)
F=A([1,3],1:2)
x=A\b

%% Question 4
I=eye(5);
B=kron(I, A)


%% Question 5
A5= normrnd(10,5,5,3)
A5(A5<10)=0;
A5(A5>=10)=1;
disp(A)

%% Question 6
M = csvread('datahw1.csv');
save('datahw1.mat','M');
load('datahw1.mat');
x1 =  M(:,3); % export
x2 = M(:,4);  % RD
y = M(:,5);  % prod
x3 = M(:,6); % capital
X=[x1 x2 x3];
% Since the matlab does not have a OLS regression function for multivariate
% regression that returns both the coefficients and their standard errors,
% I use a new function to find it. 
[b_hat, se]= ols(y, X, 1)

diary hw01.out
