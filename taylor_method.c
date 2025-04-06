#include <stdio.h>
#include <math.h>
/*IMPLEMENTATION OF THE TAYLOR METHOD FOR ODE INTEGRATION
    (Taylor method of order n)
input:
    t: time
    x: position vector
    v: velocity vector
    h: time step
    order: order of the Taylor method (number of derivatives to compute)
    taylor_coeffs: function to compute the Taylor coefficients
*/


void taylor_step(
    double *t, double *x, double *y, double h, int order,
    void (*taylor_coeffs)(double, double, double, double*, int)
) {
    double coeffs[order];
    taylor_coeffs(*t, *x, *y, coeffs, order);

    double x_new = *x;
    double y_new = *y;

    double h_pow = h;
    double prev_factor = 1.0; // 0! = 1
    for (int i = 0; i < order; i++) {
        double factor = h_pow / tgamma(i + 2);  // 1/i!
        x_new += factor * coeffs[i];
        if(i<order-1){
        y_new +=factor * coeffs[i+1]; // v = dx/dt, dv/dt = coeff[1

        } 
        h_pow *= h;
    }

    *x = x_new;
    *y = y_new;
    *t += h;
}

void map(double *x, double *y, double *t, double h,double tf, int order,char *filename,
         void (*taylor_coeffs)(double, double, double, double*, int)) {
    // Open file for gnuplot
    FILE *gnuplot_file = fopen(filename, "w");
    if (!gnuplot_file) {
        perror("Failed to open output file");
        return;
    }
    fprintf(gnuplot_file, "x\t y\n"); // header for gnuplot file
    while (*t < tf) { // Taylor method loop
        fprintf(gnuplot_file, "%.10f\t%.10f\n", x[0], x[1]); //write results
        taylor_step(t, x, y, h, order, taylor_coeffs); // perform Taylor step
        
    }
    fclose(gnuplot_file); // close file
}