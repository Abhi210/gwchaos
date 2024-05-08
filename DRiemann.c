// #include <iostream.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 4

double** cov_metric_tensor(double r, double theta, double a, double M) {
    // Define the metric tensor function.
    // Parameters:
    //     a (float): Parameter 'a' in the metric.
    //     r (float): Parameter 'r' in the metric.
    //     theta (float): Parameter 'theta' in the metric.
    // Returns:
    //     double**: The metric tensor at the given values of a, r, and theta.
    // const int N=4;
    double **g = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++)
        g[i] = (double *)malloc(N * sizeof(double));

    double Delta = a*a- 2 * M * r + r*r;
    double Sigma = r*r + a*a * cos(theta)*cos(theta) ;

    g[0][0] = -(1.0 - 2.0 * M * r / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = (
        pow(sin(theta) ,2)
        * (pow((r*r + a*a) ,2) - Delta * a*a * pow(sin(theta) ,2))
        / Sigma
    );
    g[0][3] = -2.0 * a * M * r * pow(sin(theta), 2) / Sigma;
    g[3][0] = g[0][3];

    return g;
}


double** cot_metric_tensor(double r, double theta, double a, double M) {
    // Define the metric tensor function.
    // Parameters:
    //     a (float): Parameter 'a' in the metric.
    //     r (float): Parameter 'r' in the metric.
    //     theta (float): Parameter 'theta' in the metric.
    // Returns:
    //     double**: The metric tensor at the given values of a, r, and theta.

    double **g = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++)
        g[i] = (double *)malloc(N * sizeof(double));

    double Delta = a*a- 2 * M * r + r*r;
    double Sigma = r*r + a*a * cos(theta)*cos(theta) ;

    double denom = Delta * sin(theta) * sin(theta);

    g[0][0] = (
        -(
            sin(theta) *sin(theta) 
            * ((r*r + a*a) *(r*r + a*a) - Delta * a*a * sin(theta) *sin(theta) )
            / Sigma
        )
        / denom
    );

    g[1][1] = Delta / Sigma;

    g[2][2] = 1.0 / Sigma;

    g[3][3] = (Delta - a*a * sin(theta) * sin(theta)) / (denom * Sigma);

    g[0][3] = (-2.0 * a * M * r * sin(theta) * sin(theta) / Sigma) / denom;

    g[3][0] = g[0][3];

    return g;
}

double*** kerr_christoffel(double r, double theta, double a, double M) {

    double ***cs = (double ***)malloc(N * sizeof(double **));
    for (int i = 0; i < N; i++) {
        cs[i] = (double **)malloc(N * sizeof(double *));
        for (int j = 0; j < N; j++)
            cs[i][j] = (double *)malloc(N * sizeof(double));
    }

    // Definitions
    double Delta = a*a- 2 * M * r + r*r;
    double Sigma = r*r + a*a * cos(theta)*cos(theta) ;

    cs[3][0][1] = M * (2 * r*r - Sigma) / (Delta * Sigma*Sigma) * a;

    cs[0][0][1] = (
        M * (r*r + a*a) * (r*r - a*a * cos(theta) *cos(theta)) / (Sigma*Sigma * Delta)
    );

    cs[3][0][2] = -2 * M * a * r * (cos(theta) /sin(theta)) / (Sigma*Sigma);
    cs[0][0][2] = -2 * M * a*a * r *sin(theta) * cos(theta) / (Sigma*Sigma);

    cs[0][1][3] = (
        -M
        * a
        * (2 * r*r * (r*r + a*a) + Sigma * (r*r - a*a))
        *sin(theta) *sin(theta)
        / (Delta * Sigma*Sigma)
    );

    cs[3][1][3] = (
        r * Sigma * (Sigma - 2 * M * r)
        - M * a*a * (2 * r*r - Sigma) *sin(theta) *sin(theta)
    ) / (Delta * Sigma*Sigma);

    cs[0][2][3] = M * a*a*a * r *sin(theta) *sin(theta) *sin(2 * theta) / (Sigma*Sigma);

    cs[3][2][3] = (
        (cos(theta) /sin(theta))
        * (Sigma*Sigma + 2 * M * a*a * r *sin(theta) *sin(theta))
        / (Sigma*Sigma)
    );

    cs[1][0][0] = M * Delta * (2 * r*r - Sigma) / (pow(Sigma,3));

    cs[1][0][3] = -cs[1][0][0] * a *sin(theta) *sin(theta);

    cs[2][0][0] = -M * a * r *sin(2 * theta) / (pow(Sigma,3) * a);
    cs[2][0][3] = (
        2.0 * M * r * a * (r*r + a*a) *sin(theta) * cos(theta) / (pow(Sigma,3))
    );

    cs[1][1][1] = r / Sigma + (M - r) / Delta;

    cs[1][2][2] = -r * Delta / Sigma;
    cs[2][1][2] = -cs[1][2][2] / Delta;

    cs[1][1][2]  = -(a*a) *sin(2 * theta) / (2 * Sigma);
    cs[2][2][2] = cs[1][1][2];
    cs[2][1][1] = -cs[1][1][2] / Delta;

    cs[1][3][3] = (
        -Delta
        * (r * Sigma*Sigma - M * a*a * (2 * r*r - Sigma) *sin(theta) *sin(theta))
        *sin(theta) *sin(theta)
        / pow(Sigma,3)
    );

    cs[2][3][3] = (
        -(Delta * Sigma*Sigma + 2 * M * r * (r*r + a*a) * (r*r + a*a))
        *sin(2 * theta)
        / (2 * pow(Sigma,3))
    );

    cs[0][1][0] = cs[0][0][1];
    cs[0][2][0] = cs[0][0][2];
    cs[0][3][1] = cs[0][1][3];
    cs[0][3][2] = cs[0][2][3];

    cs[1][3][0] = cs[1][0][3];
    cs[1][2][1] = cs[1][1][2];

    cs[2][3][0] = cs[2][0][3];
    cs[2][2][1] = cs[2][1][2];

    cs[3][1][0] = cs[3][0][1];
    cs[3][2][0] = cs[3][0][2];
    cs[3][3][1] = cs[3][1][3];
    cs[3][3][2] = cs[3][2][3];

    return cs;

}


int main() {
    // const int N=4;
    double r = 6.0;
    double theta = M_PI/2;
    double a = 0.8;
    double M = 1.0;

   double ***cs = kerr_christoffel(r, theta, a, M);

    printf("The Christoffel symbols cs:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                printf("%lf ", cs[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }

    // Free allocated memory
    // for (int i = 0; i < N; i++)
    //     free(g[i]);
    // free(g);

    return 0;
}


