function bootstat = BootStat1(X,W,Y,t,B,ID0,ID1,IDq1c,IDq0t,r)

     
 [N,~] = size(X);
  
  
parfor k = 1:B 
    
                S = (datasample(1:N,N))';
      bootstat(k) = BootPSM1(S,X,W,Y,t,ID0,ID1,IDq1c,IDq0t,r);
    
end;