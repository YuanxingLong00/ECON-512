function T = BootPSM1(S,X,W,Y,t,ID0,ID1,IDq1c,IDq0t,r)

% Bootstrap resample   : S
%           Covariates : X
%           Treatments : W
%           Outcomes   : Y

 [N,~] = size(W);
   N1  = sum(W);
   N0  = N-N1;
   
  
%%%Step 1: Drawing new covariates------------------------------------------
     Xs = X(S,:);                 %bootstrap covariates

 
%%%Step 2: Drawing new treatments------------------------------------------
    p1  = glmval(t,Xs,'logit','constant','off');
    Ws  = binornd(1,p1);
   

%%%Step 3: Calculating the bootstrap logit---------------------------------
    ts  = glmfit(Xs,Ws,'binomial','constant', 'off');
    ps  = glmval(ts,X,'logit','constant','off');        %Note:evaluated at X
  
  
%%%Step 4: Calculating the matching function-------------------------------
    ID0s  = knnsearch(ps(W==1), ps(W==0));  %closest match in treatment set to control values
    ID1s  = knnsearch(ps(W==0), ps(W==1));  %closest match in control set to treatment values
   
  % ID0s  = knnsearch(ps(N0+1:N), ps(1:N0));  %closest match in treatment set to control values
  % ID1s  = knnsearch(ps(1:N0), ps(N0+1:N));  %closest match in control set to treatment values
   
   %computing the matching functions
       c1 = 1:1:N1;                            
       c0 = 1:1:N0; 
      K1t = histc(ID0s,c1');                    %number of times the treatments are used as a match
      K0c = histc(ID1s,c0');                    %number of times the controls are used as a match 
  
       K0 = [K0c; K0c(IDq0t)];                %imputed matching functions
       K1 = [K1t(IDq1c); K1t];
    
   
%%%Step 5: Calculating the error terms-------------------------------------
   %(Conditional expectation using smoothing splines)
   % yhat0 = csaps(ps(1:N0),Y(1:N0),0.99,ps);  
   % yhat1 = csaps(ps(N0+1:N),Y(N0+1:N),0.99,ps);
 
   %(Conditional expectation using polynomials)
      b0 = polyfit(ps(1:N0),Y(1:N0),r);
      b1 = polyfit(ps(N0+1:N),Y(N0+1:N),r);
   yhat0 = polyval(b0, ps);
   yhat1 = polyval(b1, ps);
 
   %determining the error terms 
      e1 =  yhat1 - yhat0;                       %e_1
    e20c =  Y(1:N0) - yhat0(1:N0);               %e_2 for controls ie W=0        
    e21t =  Y(N0+1:N) - yhat1(N0+1:N);           %e_2 for treated ie W=1
 
  %1-NN matching  
     e20 = [e20c; e20c(ID1)];                   
     e21 = [e21t(ID0); e21t];
   
   %calculating the error terms using the observed samples at bootstrap ps
   %parameter
     nu1 = (1 + K1).* e21;     %constructing the net error term for original treated
     nu0 = (1 + K0).* e20;
      
   %net bootstrap errors
      es = e1(S) - (1-Ws).*nu0(S) + Ws.*nu1(S);   
   
   
%%%Step 6: Bootstrap statistic---------------------------------------------

     Xi = mean( e1 - (1-ps).*nu0 + ps.*nu1 ); %Re-centering
  
      T = mean(es)-Xi;              %bootstrap statistic
   