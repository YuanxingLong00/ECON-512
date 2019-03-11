function q = invquantile(D, p)
  % D - Data Vector
  % p - quantile probability
  % q - quantile
  [f,x]= ecdf(D);
  ind= find(f>= p, 1, 'first');
  ind1=ind-1;
  if f(ind)~= f(ind1)
      wt= (f(ind)-p)/(f(ind)-f(ind1));
      q= wt*x(ind1)+(1-wt)*x(ind);
  else
      q=x(ind);
  end
  
      