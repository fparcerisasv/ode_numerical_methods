#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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


int main(int argc, char *argv[]) { //usage: ./compute_map <x_initial> <y_initial>
    
    if (argc < 3 || argc > 4) {
        printf("Usage: %s <x_initial> <y_initial> [s]\n", argv[0]);
        printf("  s (optional): Save .dat file if provided\n");
        return 1;
    }
    //GENEREATE FLOW FROM STARTING POINT
    /*
    *
    *           REQUIRES GNUPLOT
    * 
    *  */
    //params
    double t = 0.0;
    double tf = 20.0;
    double h = 0.01;
    double tol = 1e-6;
    
    double x[2];
    x[0] = atof(argv[1]);
    x[1] = atof(argv[2]);
    int save_data = (argc == 4) ? atoi(argv[3]) : 0; // Check if 's' is provided
    int n = 2;
   

    int order = 4; // Order of the Taylor method
    

    int iterations = 100000; // total number of iterations
    int plot_iterations = 100000; // number of iterations to plot
    
    char dat_filename[100];
    char png_filename[100];
    
    snprintf(dat_filename, sizeof(dat_filename), "outputs/ode_flow_%.2f_%.2f.dat", x[0], x[1]);
    snprintf(png_filename, sizeof(png_filename), "outputs/ode_flow_%.2f_%.2f.png", x[0], x[1]);

    FILE *gnuplot_file_taylor_large = fopen(dat_filename, "w");
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
    fclose(gnuplot_file_taylor_large);
    //Plot with gnuplot
    char gnuplot_commands[6][200];
    snprintf(gnuplot_commands[0], sizeof(gnuplot_commands[0]), "set terminal pngcairo");
    snprintf(gnuplot_commands[1], sizeof(gnuplot_commands[1]), "set output \"%s\"", png_filename);
    strcpy(gnuplot_commands[2], "set title \"Map\"");
    strcpy(gnuplot_commands[3], "set xlabel \"x\"");
    strcpy(gnuplot_commands[4], "set ylabel \"y\"");
    snprintf(gnuplot_commands[5], sizeof(gnuplot_commands[5]), 
             "plot \"%s\" using 1:2 with lines title \"phi(t)\"", dat_filename);

    FILE *gnuplot_pipe = popen("gnuplot", "w");  // Removed -persistent
    for(int i = 0; i < 6; i++) {
        fprintf(gnuplot_pipe, "%s \n", gnuplot_commands[i]);
        fflush(gnuplot_pipe);  // Flush after each command
    }
    
    // Add small delay and explicit exit
    fprintf(gnuplot_pipe, "exit\n");
    fflush(gnuplot_pipe);
    
    pclose(gnuplot_pipe);
    if (!save_data) {
        remove(dat_filename); // Remove the .dat file if 's' is not provided
    }
    return 0; 
}
