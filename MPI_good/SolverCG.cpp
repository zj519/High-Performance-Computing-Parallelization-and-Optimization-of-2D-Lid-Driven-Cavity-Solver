/**
 * @file SolverCG.cpp
 * @brief Implementation of the SolverCG class for solving linear systems using the Conjugate Gradient method in parallel.
 */
#include <iostream>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cblas.h>
#include "mpi.h"
#include <omp.h>
#include "SolverCG.h"
#include "Communicator.h"

class Communicator;
using namespace std;

#define IDX(I, J) ((J) * Nx_local + (I))

/**
 * @brief Constructs a SolverCG object.
 * 
 * Initializes the SolverCG object with grid dimensions, spatial discretization,
 * and a communicator for parallel execution.
 * 
 * @param pNx Number of grid points in the x-direction.
 * @param pNy Number of grid points in the y-direction.
 * @param pdx Grid spacing in the x-direction.
 * @param pdy Grid spacing in the y-direction.
 * @param communicator Pointer to the Communicator object for MPI operations.
 */
SolverCG::SolverCG(int pNx, int pNy, double pdx, double pdy, Communicator* communicator)
    : p(communicator->p), coords(communicator->coords), communicator(communicator)
{
    dx = pdx;
    dy = pdy;
    Nx = pNx;
    Ny = pNy;
    Nx_local = communicator->Nx_local;
    Ny_local = communicator->Ny_local;
    int n = Nx*Ny;

    r = new double[n];
    k = new double[n];
    z = new double[n];
    t = new double[n]; 

    bufRecvTop = new double[Nx_local];
    bufRecvBottom = new double[Nx_local];
    bufRecvLeft = new double[Ny_local];
    bufRecvRight = new double[Ny_local];
}

/**
 * @brief Destructor for the SolverCG class.
 * 
 * Cleans up dynamically allocated memory.
 */
SolverCG::~SolverCG()
{
    delete[] r;
    delete[] k;
    delete[] z;
    delete[] t;
}

/**
 * @brief Applies the linear operator (matrix) to a vector.
 * 
 * This method simulates the multiplication of a sparse matrix, derived from the discretization
 * of a differential operator, with a vector. It accounts for boundary conditions and uses MPI
 * for parallel communication of boundary values.
 * 
 * @param in Input vector.
 * @param out Output vector, result of the operation.
 */
void SolverCG::ApplyOperator(double* in, double* out) {
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    int x_ini = communicator->x_ini;
    int end_x = communicator->end_x;
    int y_ini = communicator->y_ini;
    int end_y = communicator->end_y;
    int rank;
    MPI_Comm_rank(communicator->cart_comm, &rank);
    communicator->ExchangeBoundaryData(in, bufRecvTop, bufRecvBottom, bufRecvLeft, bufRecvRight);
    
    for (int j = y_ini; j < end_y; j ++) {
        for (int i = x_ini; i < end_x; ++i) {
            double V_left = (coords[1] > 0 && i == 0) ? bufRecvLeft[j] : in[IDX(i - 1, j)];
            double V_right = (coords[1] < p - 1 && i == Nx_local-1) ? bufRecvRight[j] : in[IDX(i + 1, j)];
            double V_bottom = (coords[0] < p - 1 && j == 0) ? bufRecvBottom[i] : in[IDX(i, j - 1)];
            double V_top = (coords[0] > 0 && j == Ny_local-1) ? bufRecvTop[i] : in[IDX(i, j + 1)];

            out[IDX(i,j)] = ( - V_left+ 2.0*in[IDX(i,   j)]-  V_right)*dx2i+ ( - V_top+ 2.0*in[IDX(i,  j)]-  V_bottom)*dy2i;
        }
    }
}

/**
 * @brief Applies the preconditioner to a vector.
 * 
 * This method simulates the application of a preconditioning step to a vector, which
 * can improve the convergence rate of the CG method. It currently uses a simple diagonal
 * preconditioner based on the inverse of the diagonal elements of the matrix.
 * 
 * @param in Input vector.
 * @param out Output vector, result of the preconditioning.
 */
void SolverCG::Precondition(double* in, double* out) {
    int i, j;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    double factor = 2.0*(dx2i + dy2i);

    int x_ini = communicator->x_ini;
    int end_x = communicator->end_x;
    int y_ini = communicator->y_ini;
    int end_y = communicator->end_y;

    for (j = y_ini; j < end_y; ++j) {
        for (i = x_ini; i < end_x; ++i) {
            out[IDX(i,j)] = in[IDX(i,j)]/factor;
        }
    }

    for (int i = 0; i < Nx_local; ++i) {
        if (coords[0] == 0)    out[IDX(i, Ny_local-1)] = in[IDX(i, Ny_local-1)];
        if (coords[0] == p - 1)  out[IDX(i, 0)] = in[IDX(i,0)];
        }
    
    for (int j = 0; j < Ny_local; ++j) {
         if (coords[1] == 0) out[IDX(0, j)] = in[IDX(0, j)];
         if (coords[1] == p - 1) out[IDX(Nx_local - 1, j)] = in[IDX(Nx_local - 1, j)];
        }
}

/**
 * @brief Imposes boundary conditions on a vector.
 * 
 * This method sets the appropriate boundary conditions for the problem being solved. It is
 * specific to the type of differential equation and the physical or geometric boundary conditions.
 * 
 * @param inout Vector on which to impose boundary conditions. Modified in place.
 */
void SolverCG::ImposeBC(double* inout) {
    // Boundaries
    for (int i = 0; i < Nx_local; ++i) {
        if (coords[0] == 0) inout[IDX(i, Ny_local-1)] = 0.0;
        if (coords[0] == p-1)  inout[IDX(i, 0)] = 0.0;
    }
    
    for (int j = 0; j < Ny_local; ++j) {
        if (coords[1] == 0) inout[IDX(0, j)] = 0.0;
        if (coords[1] == p - 1)  inout[IDX(Nx_local-1, j)] = 0.0;
    }
}

/**
 * @brief Solves the linear system Ax = b using the CG method.
 * 
 * Executes the Conjugate Gradient method to find the solution vector x that satisfies the
 * equation Ax = b, where A is a sparse matrix and b is the right-hand side vector. The method
 * uses parallelization for matrix-vector multiplication and preconditioning.
 * 
 * @param b Right-hand side vector of the linear system.
 * @param x Solution vector. Contains an initial guess on input and the solution on output.
 */
 

void SolverCG::Solve(double* b, double* x) {

    unsigned int n = Nx_local*Ny_local;
    int g;
    double alpha, alpha_glob;
    double beta, beta_glob;
    double eps;
    double tol = 0.001;

    double dot_loc = cblas_ddot(n, b, 1, b, 1);
    double dot_glob;
    MPI_Allreduce(&dot_loc, &dot_glob, 1, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
    eps = sqrt(dot_glob);

    if (eps < tol*tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps << endl;
        return;
    }
    ApplyOperator(x, t); 
    cblas_dcopy(n, b, 1, r, 1);       
    ImposeBC(r); 
    cblas_daxpy(n, -1.0, t, 1, r, 1);
    Precondition(r, z); 
    cblas_dcopy(n, z, 1, k, 1);       
    bool shouldBreak = false;
    g = 0;
        do {
            g++;
            ApplyOperator(k, t); 
            alpha = cblas_ddot(n, t, 1, k, 1);  
            MPI_Allreduce(&alpha, &alpha_glob, 1, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
            alpha = cblas_ddot(n, r, 1, z, 1) / alpha_glob; 
            MPI_Allreduce(&alpha, &alpha_glob, 1, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);

            beta = cblas_ddot(n, r, 1, z, 1);  
            MPI_Allreduce(&beta, &beta_glob, 1, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);

            cblas_daxpy(n, alpha_glob, k, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
            cblas_daxpy(n, -alpha_glob, t, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k
            dot_loc = cblas_ddot(n, r, 1, r, 1);
            MPI_Allreduce(&dot_loc, &dot_glob, 1, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
            eps = sqrt(dot_glob);
            if (eps < tol * tol) {
                    break;
                }

            Precondition(r, z);

            beta = cblas_ddot(n, r, 1, z, 1) / beta_glob;
            MPI_Allreduce(&beta, &beta_glob, 1, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
            cblas_dcopy(n, z, 1, t, 1);
            cblas_daxpy(n, beta_glob, k, 1, t, 1);
            cblas_dcopy(n, t, 1, k, 1);
            
            if(shouldBreak){ break; }
            //cout << g << endl;
        } while (g < 5000); 

        if (g == 5000) {
            cout << "FAILED TO CONVERGE" << endl;
            exit(-1);
        }
    }