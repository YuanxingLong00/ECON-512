clear
N=1000000;
se=0;
for i=1:200
f=@(x) (sqrt(1-x^2));
a=0;
b=1;
pie=4*int_NC(f,a,b,N);
se=(pie-pi)^2 +se;
i
end
mse = se/200