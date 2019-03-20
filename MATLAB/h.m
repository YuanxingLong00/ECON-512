function foc= h(x, y)

foc=1-x+x*exp(2-x)/(1+exp(2-x)+exp(2-y));
end