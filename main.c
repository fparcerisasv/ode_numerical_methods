#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RKF45.h" 

// Parameters for the ODE
const double k = 0.08;
const double b0 = 4.0;
const double b1 = 15.77;

//ODE system
void ode(double t, double *x, int n, double *f) {
    f[0] = x[1]; // dx/dt = y
    f[1] = -k * x[1] - pow(x[0], 3) + b0 + b1 * cos(t); // dy/dt
}

int main() { 
    //EXAMPLE USAGE OF RKF45

    //RKF45 parameters
    double t = 0.0;
    double tf = 20.0;
    double h = 0.01;
    double tol = 1e-6;
    double x[2] = {0.0, 0.0}; // x = 0, y = 0
    int n = 2;
    int sc = 1;
    double error;

    //open file for gnuplot
    FILE *gnuplot_file = fopen("outputs/rkf45_output.dat", "w");
    if (!gnuplot_file) {
        perror("Failed to open output file");
        return 1;
    }

    fprintf(gnuplot_file, "# t\t x\t y\t error\n");
    while (t < tf) { // RKF45 loop
        double next_tf = tf;
        int result = rkf45(&t, x, n, &h, sc, tol, &next_tf, &error, ode);
        fprintf(gnuplot_file, "%.10f\t%.10f\t%.10f\t%.10e\n", t, x[0], x[1], error); //write results
        if (result == 1) break;
    }

    fclose(gnuplot_file);
   
    //EXAMPLE USAGE OF TAYLOR METHOD

    return 0;
}
