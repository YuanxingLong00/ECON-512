function int = int_NC(f,a,b,N)
h=(b-a)/N;
X=(a+h/2):h:(b-h/2);
fv=arrayfun(f,X);
int=h* sum(fv);