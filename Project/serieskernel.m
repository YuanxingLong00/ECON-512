%% Kernel estimation for muwip0 and muwip1
K= @(u) ( 0.75*(ones(N,1)-u.^2).*indicator(u) ); % The kernel to be used 
for i=1:N
       muwip0= sum(Y.*W0.*K( (pXhattemp- pXhattemp(i))/h) )/sum(W0.*K( (pXhattemp- pXhattemp(i))/h) ); 
        muwip1= sum(Y.*W.*K( (pXhattemp- pXhattemp(i))/h) )/sum(W.*K( (pXhattemp- pXhattemp(i))/h) );
       e2Wi(i)= Y(i)-( W0(i)*muwip0+W(i)*muwip1 ) ;
        e1(i)= muwip1-muwip0-tauhats; %% use tauhats
end 