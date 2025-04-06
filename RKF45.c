#include <math.h>
#include <stdlib.h>
#include <string.h>

// VALUES FOR THE RKF45 METHOD
#define STAGES 6
double a[STAGES][STAGES] = {
    {0},
    {1.0/4.0},
    {3.0/32.0,         9.0/32.0},
    {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0},
    {439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0},
    {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}
};
double c[STAGES] = {0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0};
double b_5[STAGES] = {16.0/135.0, 0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
double b_4[STAGES] = {25.0/216.0, 0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0};

/* IMPLEMENTATION OF THE RKF45 METHOD FOR ODE INTEGRATION
   (Runge-Kutta-Fehlberg method of order 4 and 5) 
input:
    at: Time. input: time corresponding to the present initial condition.
        output: new value corresponding to the new initial condition.
    x: Position. input: present initial condition. output: new position at
        time *at.
    n: dimension of the system of odes.
    ah: time step (it can be modifed by the routine according to the given
        threshold). If negative, the integration goes backwards. On exit, it
        will contain the time step for the next rkf78 call
    sc: stepsize control.
        0: no stepsize control, the step *ah is used
        !=0: stepsize control according to the threshold tol
    tol: threshold to control the integration error
    atf: nal integration time. if NULL, it is ignored. Otherwise, if the
        stepsize is too large (*at+*ah>*atf), it is reduced so that the
        new time at is exactly atf (in that case, the function returns 1)
    aer: if NULL, the routine stops if the estimated error is larger than tol.
        if not NULL, the integration returns the estimated absolute error
        of the performed step. This allows, for instance, to integrate with
        a constant stepsize (sc= 0) and to know an estimate of the error.
    ode: pointer to the the vector eld. The parameters of the function are:
        (t,x,n,f), where t is the time, x the position vector, n the
        dimension and f the value of the vector eld at (t,x)
output:
     0: ok. 
     1: ok and at=tf
*/
int rkf45(double *at, double *x, int n, double *ah, int sc, double tol, double *atf, double *aer, void (*ode)(double,double*,int,double*)){

    double **k = malloc(STAGES * sizeof(double*));
    for (int i = 0; i < STAGES; i++) {
        k[i] = malloc(n * sizeof(double));
    }

    // Allocate memory for temporary variables and initialize them
    double *xtemp = malloc(n * sizeof(double));
    double *x4 = malloc(n * sizeof(double));
    double *x5 = malloc(n * sizeof(double));
    double *err_vec = malloc(n * sizeof(double));
    
    double t = *at;
    double h = *ah;
    int step_accepted = 0;
    int result = 0;

    while(!step_accepted){ //repeat until step is accepted
        // k1 = f(t,x) 
            ode(t, x, n, k[0]); // k1 = f(t,x)

            for (int j = 1; j < STAGES; j++) {

                // Calculate next xtemp based on previous k values
                for (int i = 0; i < n; i++) {
                    xtemp[i] = x[i];
                    for (int l= 0; l < j; l++) {
                        xtemp[i] += h * a[j][l] * k[l][i];
                    }
                }
               
                ode(t + c[j] * h, xtemp, n, k[j]); // k[j] = f(t + c[j] * h, xtemp)
            }

            //Calculate x4
            for (int i = 0; i < n; i++) {
                        double b_ = 0.0;
                for (int j = 0; j < STAGES; j++) {
                    b_ += b_4[j] * k[j][i];
                }
                x4[i] = x[i] + h * b_;
            }
            //Calculate x5
            for (int i = 0; i < n; i++) {
                double b_ = 0.0;
                for (int j = 0; j < STAGES; j++) {
                    b_ += b_5[j] * k[j][i];
                }
                x5[i] = x[i] + h * b_;
            }

            //Compute estimated error
            double err = 0.0;
            for (int i = 0; i<n; i++){
                err_vec[i] = (1.0/360.0) * k[0][i] - (128.0/4275.0) * k[2][i] - (2197.0/75240.0) * k[3][i] 
                            + (1.0/50.0) * k[4][i] + (2.0/55.0) * k[5][i];
                err += err_vec[i] * err_vec[i];
            }
            err = sqrt(err) * fabs(h); // estimated error
            if (aer) *aer = err; // store estimated error 

            double hN = 0.9 * h * pow(tol / err, 0.2); // new stepsize
            
            if(!sc ||(sc && err < tol)){ // stepsize control
               step_accepted = 1; // step accepted
               t += h; // update time
                for (int i = 0; i < n; i++) {
                     x[i] = x5[i]; // update position
                }
            } else { // error too large, reduce step size
                h = hN; // reduce step size
                if (atf && t + h > *atf) { // check if new step size is too large
                    h = *atf - t; // reduce step size to reach atf
                    result = 1; // indicate that we reached atf
                }
            }
        }

        *at = t; // update time
        *ah = h; // update step size
        
    //free allocated memory
    for (int i = 0; i < STAGES; i++) {
        free(k[i]);
    }
    free(k);
    free(xtemp);
    free(x4);
    free(x5);
    free(err_vec);
    return result;
}