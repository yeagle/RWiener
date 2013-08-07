#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double dwiener_d(double q, double alpha, double tau, double beta, double delta, int give_log)
{
  double kl, ks, ans;
  int k,K;
  double err = 1e-10;

  if (R_IsNaN(q + delta + alpha + beta + tau)) return R_NaN;
  if (!R_finite(q) || !R_finite(alpha)) return 0;
  if (beta < 0 || beta > 1 || alpha <= 0 || tau <= 0) return R_NaN;

  // q is negative for lower bound and positive for upper bound
  // extract RT and accuracy from q
  if (q<0) {
    q = fabs(q);
  }
  else {
    beta = 1-beta;
    delta = -delta;
  }
  
  q = q-tau; // remove non-decision time from q
  q = q/pow(alpha,2); // convert t to normalized time tt

  // calculate number of terms needed for large t
  if (M_PI*q*err<1) { // if error threshold is set low enough
      kl=sqrt(-2*log(M_PI*q*err)/(pow(M_PI,2)*q)); // bound
      kl=(kl>1/(M_PI*sqrt(q))) ? kl : 1/(M_PI*sqrt(q)); // ensure boundary conditions met
  }
  else { // if error threshold set too high
      kl=1/(M_PI*sqrt(q)); // set to boundary condition
  }
  // calculate number of terms needed for small t
  if ((2*sqrt(2*M_PI*q)*err)<1) { // if error threshold is set low enough
      ks=2+sqrt(-2*q*log(2*sqrt(2*M_PI*q)*err)); // bound
      ks=(ks>sqrt(q)+1) ? ks : sqrt(q)+1; // ensure boundary conditions are met
  }
  else { // if error threshold was set too high
      ks=2; // minimal kappa for that case
  }

  // compute density: f(tt|0,1,beta)
  ans=0; //initialize density
  if (ks<kl) { // if small t is better (i.e., lambda<0)
      K=ceil(ks); // round to smallest integer meeting error
      for (k=-floor((K-1)/2); k<=ceil((K-1)/2); k++) { // loop over k
          ans=ans+(beta+2*k)*exp(-(pow((beta+2*k),2))/2/q); // increment sum
      }
      ans= give_log ? log(ans)-0.5*log(2)-M_LN_SQRT_PI-1.5*log(q) : ans/sqrt(2*M_PI*pow(q,3)); // add constant term
  }
  else { // if large t is better...
      K=ceil(kl); // round to smallest integer meeting error
      for (k=1; k<=K; k++) {
          ans=ans+k*exp(-(pow(k,2))*(pow(M_PI,2))*q/2)*sin(k*M_PI*beta); // increment sum
      }
      ans= give_log ? log(ans)+2*M_LN_SQRT_PI : ans*M_PI; // add constant term
  }

  // convert to f(t|v,a,w) and return result
  return give_log ? 
    ans+((-delta*alpha*beta -(pow(delta,2))*(q*pow(alpha,2))/2)-log(pow(alpha,2))) : 
    ans*exp(-delta*alpha*beta -(pow(delta,2))*(q*pow(alpha,2))/2)/(pow(alpha,2));
}

SEXP dwiener(SEXP q, SEXP alpha, SEXP tau, SEXP beta, SEXP delta, SEXP give_log) {
  double d;
  SEXP value;

  if (fabs(REAL(q)[0]) <= REAL(tau)[0]) {
    if(LOGICAL(give_log)[0]) d = -1.0/0.0; // -inf
    else d = 0;
  }
  else {
    d = dwiener_d(REAL(q)[0], REAL(alpha)[0], REAL(tau)[0], REAL(beta)[0], REAL(delta)[0], LOGICAL(give_log)[0]);
  }

  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}
