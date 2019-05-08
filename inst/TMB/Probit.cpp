#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //data
  DATA_MATRIX(X);
  DATA_VECTOR(y); 
 
  //parameters
  PARAMETER_VECTOR(para);
  

  using namespace density;

  Type nll=0.0;     // Negative log likelihood function

  int n = X.rows();
  
 
   vector<Type> Xbeta = X * para;
 
   
  for(int i=0;i<n;i++)
  {
     Type eta = Xbeta(i);
     Type prob =  pnorm(eta, Type(0), Type(1));
     prob = squeeze(prob);
     nll -= dbinom(y(i), Type(1), prob, true);
   }
   return nll;
}
