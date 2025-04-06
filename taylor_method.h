#ifndef TAYLOR_METHOD_H
#define TAYLOR_METHOD_H

// IMPLEMENTATION OF THE TAYLOR METHOD FOR ODE INTEGRATION
void taylor_step(
    double *t, double *x, double *y, double h, int order,
    void (*taylor_coeffs)(double, double, double, double*, int)
);
void map(double *x, double *y, double *t, double h,double tf, int order,char *filename,

    void (*taylor_coeffs)(double, double, double, double*, int));
#endif // TAYLOR_METHOD_H
