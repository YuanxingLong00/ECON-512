function bootstat = PSM_Bootstrap(X,W,Y,Q,r,B)

% Inputs: X-Covariates (Nxk) (Note: constants have to be added in X)
%         W-Treatments (Nx1)
%         Y-Outcomes   (Nx1)  
%         Q-Number of quantile partitions (Q has to be greater than 1; recommendation is Q=4, equivalent to 5 blocks)
%         r-Number of polynomials for estimating errors (recommendation is 1)
%         B-Number of bootstrap replications 
%Output:  bootstat-vector of bootstrap statistics
%Notes: 1) As of the moment the code only applies for single matches
     %  2) This code only does bootstrap for a single imputation (for now)
     %  3) This function requires the additional m files: BootStat1.m and
     %     BootPSM1.m


%-------------------Step 1: Estimate Pouplation quantities-----------------
        t  =  glmfit(X,W,'binomial','constant', 'off');   %estimated propensity score
     [N,~] = size(X);
     
     %Logit estimates of the propensity score parameter
      phat =  glmval(t,X,'logit','constant','off'); 
   
     
%-------------------Step 2: Imputation for the `within' errors-------------
         p = mean(W==1);     
       Cov = p*nancov(X(W==1,:)) + (1-p).*nancov(X(W==0,:));
          
       ID0 = knnsearch(X(W==1,:), X(W==0,:), 'Distance', 'mahalanobis', 'Cov', Cov);  %closest match in treatment set to control values
       ID1 = knnsearch(X(W==0,:), X(W==1,:), 'Distance', 'mahalanobis', 'Cov', Cov);  %closest match in control set to treatment values

    
%-------------------Step 3: Stratify observations in blocks----------------
  
        q = quantile(phat, Q);           %careful when choosing Q=1!
    [~,L] = histc(phat, [0,q,1]);   
       L0 = L(W==0);
       L1 = L(W==1);
 
      N1  = sum(W);
      N0  = N-N1;
   
       c1 = 1:1:N1;                            
       c0 = 1:1:N0;
    
    IDq0t = zeros(N1,1);
    IDq1c = zeros(N0,1);
        
    for l = 1:Q+1
     
        Nl1 = sum(L1 == l);              %number of treated terms in stratum l
        N0l = sum(L0 == l);
     
       if and(Nl1~=0, N0l~=0)   
          IDq0t(L1 == l) = datasample( c0(L0 == l) , Nl1);
          IDq1c(L0 == l) = datasample( c1(L1 == l) , N0l);     
         
       elseif Nl1 == 0  
          IDq1c(L0 == l) = datasample( c1 , N0l);
                 
       else  
          IDq0t(L1 == l) = datasample( c0, Nl1);                        
       end   
      
    end    
    
   
%----------------------Step 4: Apply bootstrap algorithm-------------------   
       
   bootstat = BootStat1(X,W,Y,t,B,ID0,ID1,IDq1c,IDq0t,r);   