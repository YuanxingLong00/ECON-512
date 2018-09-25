function J =jacob(p)
x=p(1);
y=p(2);
J=[ exp(2 - x)/(exp(2 - x) + exp(2 - y) + 1) - (x*exp(2 - x))/(exp(2 - x) + exp(2 - y) + 1) + (x*exp(4 - 2*x))/(exp(2 - x) + exp(2 - y) + 1)^2 - 1,                                                                                      (x*exp(2 - x)*exp(2 - y))/(exp(2 - x) + exp(2 - y) + 1)^2; 
                                                                                      (y*exp(2 - x)*exp(2 - y))/(exp(2 - x) + exp(2 - y) + 1)^2, exp(2 - y)/(exp(2 - x) + exp(2 - y) + 1) - (y*exp(2 - y))/(exp(2 - x) + exp(2 - y) + 1) + (y*exp(4 - 2*y))/(exp(2 - x) + exp(2 - y) + 1)^2 - 1];
end
