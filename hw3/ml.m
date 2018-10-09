% Define the MLE objective function.
function mlobj=ml(beta)
load('hw3.mat');
temp= X*[beta(1);beta(2);beta(3);beta(4);beta(5);beta(6)];
temp1= -exp(temp)+y.* temp;  % drop the constant 
mlobj=  -sum(temp1, 1);
