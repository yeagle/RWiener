#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double fl01(double t, double beta, double lambda, double kappa) {
  double res = 0;

  if(lambda < 0) { // ST = small times
    for (int k=-ceil((kappa-1)/2); k<=floor((kappa-1)/2); k++) {
      res += ( (beta + 2*k)*exp(-(1.0L/2.0L)*pow(beta + 2*k, 2)/t) );
    }
    res = (1.0L/2.0L)*sqrt(2)*res/(sqrt(M_PI)*sqrt(pow(t, 3)));
  }
  else { // LT = large times
    for (int k=1; k<=ceil(kappa); k++) {
      res +=  ( k*exp(-(1.0L/2.0L)*pow(M_PI, 2)*pow(k, 2)*t)*sin(M_PI*beta*k) );
    }
    res = M_PI*res;
  }

  return res;
}

double scl01tau(double t, double beta, double lambda, double kappa) {
  double res = 0;
  double res2 = 0;

  if (lambda < 0) { // ST = small times
    for (int k=-ceil((kappa-1)/2); k<=floor((kappa-1)/2); k++) {
      res += ( -2*pow(beta + 2*k, 3)*exp(-pow(beta + 2*k, 2)
        / (2*t))/pow(2*t, 2) );
      res2 += ( (beta + 2*k)*exp(-pow(beta + 2*k, 2)/(2*t)) );
    }
    res2 = (3.0L/4.0L)*sqrt(2)*res2/(sqrt(M_PI)*(t)*sqrt(pow(t, 3)));
    res = res2 + (1.0L/2.0L)*sqrt(2)*res/(sqrt(M_PI)*sqrt(pow(t, 3)));
  }
  else { // LT = large times
    for (int k=1; k<=ceil(kappa); k++) {
      res += ( (1.0L/2.0L)*pow(M_PI, 2)*pow(k, 3)*exp(-(1.0L/2.0L)
        *pow(M_PI, 2)*pow(k, 2)*t)*sin(M_PI*beta*k) );
    }
    res = res*M_PI;
  }

  return res;
}

double scl01beta(double t, double beta, double lambda, double kappa) {
  double res = 0;

  if(lambda < 0) { // ST = small times
    for (int k=-ceil((kappa-1)/2); k<=floor((kappa-1)/2); k++) {
      res += ( exp(-(1.0L/2.0L)*pow(beta + 2*k, 2)/t) 
        - (1.0L/2.0L)*(beta + 2*k)*(2*beta + 4*k)
        *exp(-(1.0L/2.0L)*pow(beta + 2*k, 2)/t)/t );
    }
    res = (1.0L/2.0L)*sqrt(2)*res/(sqrt(M_PI)*sqrt(pow(t, 3)));
  }
  else { // LT = large times
    for (int k=1; k<=ceil(kappa); k++) {
      res += ( M_PI*pow(k, 2)*exp(-(1.0L/2.0L)*pow(M_PI, 2)*pow(k, 2)*t)*cos(M_PI*beta*k) );
    }
    res = res*M_PI;
  }

  return res;
}

SEXP kappaLT(SEXP t) {
  SEXP value;
  double err = 1e-10;

  PROTECT(value = allocVector(REALSXP, 1));

  REAL(value)[0] = ( sqrt(2)*sqrt(-log(M_PI*err*REAL(t)[0])/REAL(t)[0])/M_PI );

  UNPROTECT(1);

  return value;
}

SEXP kappaST(SEXP t) {
  double err = 1e-10;
  SEXP value;

  PROTECT(value = allocVector(REALSXP, 1));

  REAL(value)[0] = ( sqrt(2)*sqrt(-REAL(t)[0]*log(2*sqrt(2)*sqrt(M_PI)*err*sqrt(REAL(t)[0]))) + 2 );

  UNPROTECT(1);

  return value;
}

SEXP sclalpha(SEXP t, SEXP alpha, SEXP tau, SEXP beta, SEXP delta, SEXP lambda, SEXP kappa) {
  double err = 1e-10;
  SEXP value;
  double tx = REAL(t)[0] - REAL(tau)[0];

  PROTECT(value = allocVector(REALSXP, 1));

  REAL(value)[0] =  ( -(REAL(beta)[0])*(REAL(delta)[0])
    * fl01(tx/pow((REAL(alpha)[0]),2), (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0])) 
    * exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0]) 
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)/pow((REAL(alpha)[0]), 2) 
    - 2*fl01(tx/pow((REAL(alpha)[0]),2), (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0]))
    * exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0]) 
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)/pow((REAL(alpha)[0]), 3) 
    - 0 // 2*tx*exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0])-(1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx) * 0 
    / pow((REAL(alpha)[0]), 5) );

  UNPROTECT(1);

  return value;
}

SEXP scltau(SEXP t, SEXP alpha, SEXP tau, SEXP beta, SEXP delta, SEXP lambda, SEXP kappa) {
  double err = 1e-10;
  SEXP value;
  double tx = REAL(t)[0] - REAL(tau)[0];

  PROTECT(value = allocVector(REALSXP, 1));

  REAL(value)[0] = ( (1.0L/2.0L)*pow((REAL(delta)[0]), 2)
    * fl01(tx/pow((REAL(alpha)[0]),2), (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0]))
    * exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0]) 
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)/pow((REAL(alpha)[0]), 2) 
    - exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0]) 
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)
    * scl01tau((tx/pow((REAL(alpha)[0]),2)), (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0]))
    / pow((REAL(alpha)[0]), 4) );

  UNPROTECT(1);

  return value;
}

SEXP sclbeta(SEXP t, SEXP alpha, SEXP tau, SEXP beta, SEXP delta, SEXP lambda, SEXP kappa) {
  double err = 1e-10;
  SEXP value;
  double tx = REAL(t)[0] - REAL(tau)[0];

  PROTECT(value = allocVector(REALSXP, 1));

  REAL(value)[0] = ( -(REAL(delta)[0])
    * fl01(tx/pow((REAL(alpha)[0]),2), (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0]))
    * exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0]) 
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)/(REAL(alpha)[0]) 
    + exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0]) 
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)
    * scl01beta(tx, (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0]))
    / pow((REAL(alpha)[0]), 2) );

  UNPROTECT(1);

  return value;
}

SEXP scldelta(SEXP t, SEXP alpha, SEXP tau, SEXP beta, SEXP delta, SEXP lambda, SEXP kappa) {
  double err = 1e-10;
  SEXP value;
  double tx = REAL(t)[0] - REAL(tau)[0];

  PROTECT(value = allocVector(REALSXP, 1));

  REAL(value)[0] = ( (-(REAL(alpha)[0])*(REAL(beta)[0]) 
    - (REAL(delta)[0])*tx)
    * fl01(tx/pow((REAL(alpha)[0]),2), (REAL(beta)[0]), (REAL(lambda)[0]), (REAL(kappa)[0]))
    * exp(-(REAL(alpha)[0])*(REAL(beta)[0])*(REAL(delta)[0])
    - (1.0L/2.0L)*pow((REAL(delta)[0]), 2)*tx)
    / pow((REAL(alpha)[0]), 2) );

  UNPROTECT(1);

  return value;
}

