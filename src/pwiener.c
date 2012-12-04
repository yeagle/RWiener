#include <R.h>
#include <Rmath.h>

double prob_upperbound(double v, double a, double w)
{
  double e = exp(-2.0 * v * a * (1.0-w));
  if(e == R_PosInf) return 1;
  else if(v == 0 || w == 1) return (1-w);
  else return ((1 - e) / (exp(2.0*v*a*w) - e));
}

double exp_pnorm(double a, double b)
{
  double r;
  if (R_IsNaN(r) && b < -5.5) r = 1/sqrt(2) * exp(a - b*b/2) * (0.5641882/b/b/b - 1/b/sqrt(M_PI));
  else r = exp(a) * pnorm(b,0,1,1,0);
  return r;
}

int K_large(double q, double v, double a, double w)
{
  double err = 1e-10;
  double sqrtL1 = sqrt(1/q) * a/M_PI;
  double sqrtL2 = sqrt(fmax(1.0, -2/q*a*a/M_PI/M_PI * (log(err*M_PI*q/2 * (v*v + M_PI*M_PI/a/a)) + v*a*w + v*v*q/2)));
  return ceil(fmax(sqrtL1, sqrtL2));
}

int K_small(double q, double v, double a, double w, double epsilon)
{
  if(v == 0) return ceil(fmax(0.0, w/2 - sqrt(q)/2/a * qnorm(fmax(0.0, fmin(1.0, epsilon/(2-2*w))),0,1,1,0)));
  if(v > 0) return(K_small(q, -v, a, w, exp(-2*a*w*v)*epsilon));
  double S2 = w - 1 + 0.5/v/a * log(epsilon/2 * (1-exp(2*v*a)));
  double S3 = (0.535 * sqrt(2*q) + v*q + a*w)/2/a;
  double S4 = w/2 - sqrt(q)/2/a * qnorm(fmax(0.0, fmin(1.0, epsilon * a / 0.3 / sqrt(2*M_PI*q) * exp(v*v*q/2 + v*a*w))),0,1,1,0);
  return ceil(fmax(fmax(fmax(S2, S3), S4), 0.0));
}

double Fl_lower(double q, double v, double a, double w, int K)
{
  double F=0;
  for(int k=K; k>=1; k--) F = F - k / (v*v*1.0 + k*k*M_PI*M_PI/(a*1.0)/a) * exp(-v*a*w*1.0 - 0.5*v*v*q - 0.5*k*k*M_PI*M_PI/(a*1.0)/a*q) * sin(M_PI*k*w);
  return prob_upperbound(v, a, w) + 2.0*M_PI/(a*1.0)/a * F;
}

double Fs0_lower(double q, double a, double w, int K)
{
  double F=0;
  for(int k=K; k>=0; k--) {
    F = F - pnorm((-2*k - 2 + w)*a/sqrt(q),0,1,1,0) + pnorm((-2*k - w)*a/sqrt(q),0,1,1,0);
  }

  return 2*F;
}

double Fs_lower(double q, double v, double a, double w, int K)
{
  if (v == 0) return(Fs0_lower(q, a, w, K));
  double S1=0,S2=0;
  double sqt = sqrt(q);
  for(int k=K; k>=1; k--) {
    S1 = S1 + exp_pnorm(2*v*a*k, -sign(v)*(2*a*k+a*w+v*q)/sqt) -
           exp_pnorm(-2*v*a*k - 2*v*a*w, sign(v)*(2*a*k+a*w-v*q)/sqt);
    S2 = S2 + exp_pnorm(-2*v*a*k, sign(v)*(2*a*k-a*w-v*q)/sqt) -
           exp_pnorm(2*v*a*k - 2*v*a*w, -sign(v)*(2*a*k-a*w+v*q)/sqt);
  }
  return prob_upperbound(v, a, w) + sign(v) * ((pnorm(-sign(v) * (a*w+v*q)/sqt,0,1,1,0) -
           exp_pnorm(-2*v*a*w, sign(v) * (a*w-v*q)/sqt)) + S1 + S2);
}

double F_lower(double q, double v, double a, double w)
{
  /*  double sigma = 1;
      a = a / sigma;
      v = v / sigma; */
  double err = 1e-10;
  double F;
  int K_l = K_large(q, v, a, w);
  int K_s = K_small(q, v, a, w, err);
  if (K_l < 10*K_s) F = Fl_lower(q, v, a, w, K_l);
  else F = Fs_lower(q, v, a, w, K_s);
  return F;
}

double pwiener(double q, double alpha, double tau, double beta, double delta)
{
  double p;

  if(!R_finite(q)) return R_PosInf;
  if (R_IsNaN(q)) return R_NaN;
  if (fabs(q) <= tau) return 0;

  if (q < 0) { // lower boundary 0
    p = F_lower(fabs(q)-tau, delta, alpha, beta);
  }
  else { // upper boundary a
    p = F_lower(q-tau, (-delta), alpha, (1-beta));
  }

  return p;
}


double pwiener_full(double q, double alpha, double tau, double beta, double delta)
{
  if (q < 0) return R_NaN;
  if(!R_finite(q)) return R_PosInf; // infinity
  return (pwiener(q, alpha,tau,beta,delta) + pwiener(-q, alpha,tau,beta,delta));
}

