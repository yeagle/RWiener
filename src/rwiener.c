#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double r_random_walk(double alpha, double tau, double beta, double delta)
{
  double dt=0.0001;
  double t,sigma=1;
  double p = .5 * (1+((delta*sqrt(dt))/sigma));
  double a;
  //double q = .5 * (1-((mu*sqrt(dt))/sigma));
  int i = 0;
  double y = beta*alpha;

  while(y < alpha && y > 0)
  {
    GetRNGstate();
    a = unif_rand();
    PutRNGstate();
    if(a <= p) y = y + sigma*sqrt(dt);
    else y = y - sigma*sqrt(dt);
    i++;
  }
  if(y >= alpha) t = (i*dt+tau);
  else t = -(i*dt+tau);

  return t;
}

double r_rejection_based(double a, double ter, double z, double v)
{
  /*  mere copy of wdmrnd.cpp by JV, only changes:
   *  - return value double instead of void
   *  - removed *t and *x, instead returning t or -t 
   *  - added variable t (double)
   *  - replaced GNU gsl with unif_rand()
   *  - absol replaced with fabs
   *  - amin replaced with fmin
   *  - pi replaced with M_PI
   */
  double dt=1e-15,tau=.1,D=.005,totaltime,startpos,ndrt,
  zz,Aupper,Alower,radius,lambda,F,prob,tt,dir_,l,s1,s2,tnew,t_delta;
  int uu,i;
  int finish;
  double t, r;
  
  a/=10;
  z/=10;
  v/=10;

  finish = 0;
  totaltime=0;
  startpos=0;
  Aupper=a-z;
  Alower=-z;
  radius=fmin(fabs(Aupper),fabs(Alower));
  
  while (!finish) {
    if (v==0){
      lambda = 0.25*D*M_PI*M_PI/(radius*radius);
      F=1;
      prob = .5;
    } 
    else {
      lambda = 0.25*v*v/D + 0.25*D*M_PI*M_PI/(radius*radius);
      F=D*M_PI/(radius*v);
      F=F*F/(1+F*F);
      prob=exp(radius*v/D);
      prob=prob/(1+prob);
    }
    GetRNGstate();
    r = unif_rand();
    PutRNGstate();
    dir_= r<prob ? 1 : -1;
    l=-1;
    s2=0;
    
    while (s2>l) {
      GetRNGstate();
      s2 = unif_rand();
      PutRNGstate();
      GetRNGstate();
      s1 = unif_rand();
      PutRNGstate();
      tnew=0;
      t_delta=0;
      uu=0;
      
      while ( (fabs(t_delta)>dt) | (!uu) ) {
        tt = 2*++uu+1;
        t_delta = tt * (uu%2?-1:1) * pow(s1,(F*tt*tt));
        tnew += t_delta;
      }
      
      l = 1 + pow(s1,-F) * tnew;
    }/*end while (s2>l) */
    
    totaltime+=fabs(log(s1))/lambda;
    dir_=startpos+dir_*radius;
    
    if (dir_+dt>Aupper) {
      //*t=totaltime+ter;
      //*x=1;
      t = totaltime+ter;
      finish=1;
      return t;
    }
    else {
      if (dir_-dt<Alower) {
        //*t=totaltime+ter;
        //*x=0;
        t = -(totaltime+ter);
        finish=1;
        return t; 
      }
      else {
        startpos=dir_;
        radius=fmin(fabs(Aupper-startpos),fabs(Alower-startpos));
      }
    }
  } /*end while (!finish) */
}

double rwiener_d(double alpha, double tau, double beta, double delta)
{
  return r_rejection_based(alpha, tau, beta*alpha, delta);
}

SEXP rwiener(SEXP alpha, SEXP tau, SEXP beta, SEXP delta) {
  double r;
  SEXP value;

  r =  rwiener_d(REAL(alpha)[0], REAL(tau)[0], REAL(beta)[0], REAL(delta)[0]);

  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = r;
  UNPROTECT(1);
  return value;
}
