#ifndef RKF45_H
#define RKF45_H

// RKF45 integrator
int rkf45(double *at, double *x, int n, double *ah, int sc, double tol,
          double *atf, double *aer,
          void (*ode)(double, double*, int, double*));

#endif // RKF45_H
