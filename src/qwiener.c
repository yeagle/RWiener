#include <R.h>
#include <Rmath.h>

// prototypes
double pwiener(double q, double alpha, double tau, double beta, double delta);
double pwiener_full(double q, double alpha, double tau, double beta, double delta);

double qwiener_full(double p, double alpha, double tau, double beta, double delta)
{
  double pmin,pmax,pmid;
  double qmin,qmax,q;

  if (p > 1) return R_NaN;

  q = 1;
  qmin = 0;
  pmin = 0;
  qmax = R_PosInf;
  pmax = 1;

  int c=0;
  do {
    c++;
    pmid = pwiener_full(q, alpha,tau,beta,delta);
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
  } while(fabs(p-pmid) > 1e-10 && c < 1000); // defines the accuracy

  return q;
}

double qwiener(double p, double alpha, double tau, double beta, double delta)
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
    if (p>=0) pmid = pwiener(q, alpha,tau,beta,delta);
    else pmid = pwiener(-q, alpha,tau,beta,delta);
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
