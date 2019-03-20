% Define the MLE objective function.
function mlobj=ml(beta)
load('hw3.mat');
temp= X*beta;
temp1= -exp(temp)+y.* temp;  % drop the constant 
mlobj= sum(temp1, 1);
