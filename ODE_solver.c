#include "ODE_solver.h"

void RK45_stepping(double* t, double*v, double* dt, 
    ODE_func f, ODE_parameters* p, int n) {
    /*
    Take one RK4 step. Return updated solution and time.
    f: Right-hand-side function: dv/dt = f(v)
    v: current solution
    t: current time
    dt: time step
    n: number of equations
    */
    double* k1 = (double*)malloc(n * sizeof(double));
    double* k2 = (double*)malloc(n * sizeof(double));
    double* k3 = (double*)malloc(n * sizeof(double));
    double* k4 = (double*)malloc(n * sizeof(double));
    double* v_temp = (double*)malloc(n * sizeof(double));
    if (!k1 || !k2 || !k3 || !k4 || !v_temp) {
        fprintf(stderr, "Memory allocation faild\n");
        exit(EXIT_FAILURE);
    }

    // --- k1 = dt * f(v) ---
    // copy current solution to p->v as input for f
    cblas_dcopy(n, v, 1, p->v, 1);
    // call f and then the result is stored in p->v
    f((void*) p);
    // Copt the derivative into k1 and scale by dt
    cblas_dcopy(n, p->v, 1, k1, 1);
    cblas_dscal(n, *dt, k1, 1);

    // --- k2 = dt * f(v + 0.5 * k) ---
    // v_temp = v + 0.5 * k1
    cblas_dcopy(n, v, 1, v_temp, 1);
    // cblas_daxpy(n, alpha, X, incX, Y, incY): 
    // Computes Y = alpha * X + Y for an n-element vector.
    cblas_daxpy(n, 0.5, k1, 1, v_temp, 1);
    // copy v_temp to p->v as input for f
    cblas_dcopy(n, v_temp, 1, p->v, 1);
    // call f and then the result is stored in p->v
    f((void*) p);
    // Copt the derivative into k2 and scale by dt
    cblas_dcopy(n, p->v, 1, k2, 1);
    cblas_dscal(n, *dt, k2, 1);


    // Compute k3 = dt * f(v + 0.5 * k2)
    // tmp = v + 0.5 * k2
    cblas_dcopy(n, v, 1, v_temp, 1);
    cblas_daxpy(n, 0.5, k2, 1, v_temp, 1);
    // copy v_temp to p->v as input for f
    cblas_dcopy(n, v_temp, 1, p->v, 1);
    // call f and then the result is stored in p->v
    f((void*) p);
    // Copt the derivative into k3 and scale by dt
    cblas_dcopy(n, p->v, 1, k3, 1);
    cblas_dscal(n, *dt, k3, 1);

    // Compute k4 = dt * f(v + k3)
    // tmp = v + k3
    cblas_dcopy(n, v, 1, v_temp, 1);
    cblas_daxpy(n, 1.0, k3, 1, v_temp, 1);
    // copy v_temp to p->v as input for f
    cblas_dcopy(n, v_temp, 1, p->v, 1);
    // call f and then the result is stored in p->v
    f((void*) p);
    // Copt the derivative into k4 and scale by dt
    cblas_dcopy(n, p->v, 1, k4, 1);
    cblas_dscal(n, *dt, k4, 1);

    // Update solution
    // v = v + (k1 + 2*k2 + 2*k3 + k4) / 6
    // First, compute tmp = k1 + 2*k2 + 2*k3 + k4
    cblas_dcopy(n, v, 1, v_temp, 1);
    cblas_daxpy(n, 2.0, k2, 1, v_temp, 1);
    cblas_daxpy(n, 2.0, k3, 1, v_temp, 1);
    cblas_daxpy(n, 1.0, k4, 1, v_temp, 1);
    // Scale by 1/6
    cblas_dscal(n, 1.0/6.0, v_temp, 1);
    // ADD v_temp to v to update the solution
    cblas_daxpy(n, 1.0, v_temp, 1, v, 1);

    //free memory
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(v_temp);
}