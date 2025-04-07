#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RKF45.h" 
#include "taylor_method.h"
#define PI 3.14159265358979323846
// Parameters for the ODE
const double k = 0.08;
const double b0 = 4.0;
const double b1 = 15.77;

//ODE system
void ode(double t, double *x, int n, double *f) {
    f[0] = x[1]; // dx/dt = y
    f[1] = -k * x[1] - pow(x[0], 3) + b0 + b1 * cos(t); // dy/dt
}

void taylor_coeffs(double t, double x, double y, double *coeffs, int order) { //taylor coefficients for example ODE
    if (order >= 1)
    coeffs[0] = y; 

    if (order >= 2)
        coeffs[1] = -k * y + b0 - x * x * x + b1 * cos(t);

    if (order >= 3)
        coeffs[2] = -k * coeffs[1] - 3 * x * x * y - b1 * sin(t);

    if (order >= 4)
        coeffs[3] = -k * coeffs[2] - 3 * (2 * x * y * y + x * x * coeffs[1]) - b1 * cos(t);
}


int main() { 
    //1. EXAMPLE USAGE OF RKF45

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
   
    //2. EXAMPLE USAGE OF TAYLOR METHOD
    int order = 4; // Order of the Taylor method
    FILE *gnuplot_file_taylor = fopen("outputs/taylor_output.dat", "w");
    if (!gnuplot_file_taylor) {
        perror("Failed to open output file");
        return 1;
    }
    // reset variables for Taylor method
    t = 0.0; 
    tf = 20.0; 
    h = 0.01; 
    x[0] = 0.0; // x = 0
    x[1] = 0.0; // y = 0

    fprintf(gnuplot_file_taylor, "# t\t x\t y\n");
    while (t < tf) { // Taylor method loop
        fprintf(gnuplot_file_taylor, "%.10f\t%.10f\t%.10f\n", t, x[0], x[1]); //write results
        taylor_step(&t, &x[0], &x[1], h, order, taylor_coeffs);
        
    }
    fclose(gnuplot_file_taylor);

    // 3. EXAMPLE USAGE OF TAYLOR METHOD WITH MAP FUNCTION
    double x_map[2] = {0.0,0.0}; // initial x position
    double t_map = 0.0; // initial time
    double h_map = 0.01; // time step
    double tf_map = 2*PI; // final time
    char *filename = "outputs/taylor_map_output.dat"; // output file name
    int order_map = 4; // Order of the Taylor method
    
    map(&x_map[0], &x_map[1], &t_map, h_map, tf_map, order_map, filename, taylor_coeffs); // call map function
    
    // 4. RUNNING 100100 ITERATIONS AND PLOTTING LAST 10000 USING taylor method

    int iterations = 100100; // total number of iterations
    int plot_iterations = 10000; // number of iterations to plot

    // reset variables for Taylor method
    t = 0.0; 
    tf = iterations * h; // final time
    h = 0.01; 
    x[0] = 0.0; // x = 0
    x[1] = 0.0; // y = 0
    FILE *gnuplot_file_taylor_large = fopen("outputs/taylor_large_output.dat", "w");
    if (!gnuplot_file_taylor_large) {
        perror("Failed to open output file");
        return 1;
    }
    fprintf(gnuplot_file_taylor_large, "#x\t y\n");
    while (t < tf) { // Taylor method loop
        if (t > tf - plot_iterations * h) { // only write last 10000 iterations
            fprintf(gnuplot_file_taylor_large, "%.10f\t%.10f\n", x[0], x[1]); //write results
        }
        taylor_step(&t, &x[0], &x[1], h, order, taylor_coeffs);
    }
    return 0; 
}
