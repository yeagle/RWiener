#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

// prototypes
double pwiener_d(double q, double alpha, double tau, double beta, double delta);
double pwiener_full_d(double q, double alpha, double tau, double beta, double delta);

double qwiener_full_d(double p, double alpha, double tau, double beta, double delta)
{
  double pmid;
  double qmin,qmax,q;

  if (p > 1) return R_NaN;

  q = 1;
  qmin = 0;
  qmax = R_PosInf;

  int c=0;
  do {
    c++;
    pmid = pwiener_full_d(q, alpha,tau,beta,delta);
    if (fabs(p)<=pmid) { // near lower point
      qmax = q;
      q = qmin + (qmax-qmin)/2;
    }
    else { // near upper point
      qmin = q;
      if (R_finite(qmax)) 
        q = qmin + (qmax-qmin)/2;
      else
        q = q*10;
    }
    if(R_IsNaN(pmid)) return R_NaN;
    if(q>=1e+10) return R_PosInf;
  } while(fabs(p-pmid) > 1e-10 && c < 1000); // defines the accuracy

  return q;
}

double qwiener_d(double p, double alpha, double tau, double beta, double delta)
{
  double pmin,pmax,pmid;
  double qmin,qmax,q;

  if (fabs(p) > 1) return R_NaN;

  q = 1;
  qmin = 0;
  pmin = 0;
  qmax = R_PosInf;
  pmax = 1;

  int c=0;
  do {
    c++;
    if (p>=0) pmid = pwiener_d(q, alpha,tau,beta,delta);
    else pmid = pwiener_d(-q, alpha,tau,beta,delta);
    if (fabs(p)<=pmid) { // near lower point
      pmax = pmid;
      qmax = q;
      q = qmin + (qmax-qmin)/2;
    }
    else { // near upper point
      pmin = pmid;
      qmin = q;
      if (R_finite(qmax)) 
        q = qmin + (qmax-qmin)/2;
      else
        q = q*10;
    }
    if(R_IsNaN(pmid)) return R_NaN;
    if(q>=1e+10) return R_PosInf;
  } while(fabs(fabs(p)-pmid) > 1e-10 && c < 1000); // defines the accuracy

  return q;
}

SEXP qwiener(SEXP p, SEXP alpha, SEXP tau, SEXP beta, SEXP delta) {
  double q;
  SEXP value;

  q =  qwiener_d(REAL(p)[0], REAL(alpha)[0], REAL(tau)[0], REAL(beta)[0], REAL(delta)[0]);

  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = q;
  UNPROTECT(1);
  return value;
}

SEXP qwiener_full(SEXP p, SEXP alpha, SEXP tau, SEXP beta, SEXP delta) {
  double q;
  SEXP value;

  q = qwiener_full_d(REAL(p)[0], REAL(alpha)[0], REAL(tau)[0], REAL(beta)[0], REAL(delta)[0]);

  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = q;
  UNPROTECT(1);
  return value;
}
