/**
 * @file SolverCG.h
 * @brief Defines the SolverCG class for solving linear systems using the Conjugate Gradient method.
 */

#pragma once
#include "Communicator.h"

/**
 * @class SolverCG
 * @brief Implements the Conjugate Gradient method for solving linear systems in a distributed manner.
 *
 * SolverCG is designed to solve sparse linear systems that arise from discretizing
 * partial differential equations. It works in parallel by distributing the computation
 * across MPI processes using a Communicator object.
 */
class SolverCG
{
public:
    /**
     * @brief Constructs a SolverCG object with the specified grid dimensions and communicator.
     * 
     * @param pNx Number of grid points in the x-direction.
     * @param pNy Number of grid points in the y-direction.
     * @param pdx Grid spacing in the x-direction.
     * @param pdy Grid spacing in the y-direction.
     * @param communicator Pointer to the communicator for parallel operations.
     */
    SolverCG(int pNx, int pNy, double pdx, double pdy, Communicator* communicator);

    /**
     * @brief Destroys the SolverCG object.
     * 
     * Cleans up dynamically allocated memory used during the solution process.
     */
    ~SolverCG();

    /**
     * @brief Solves the linear system Ax = b using the Conjugate Gradient method.
     * 
     * @param b The right-hand side vector of the linear system.
     * @param x The solution vector to be computed. Should be initialized to an initial guess.
     */
    void Solve(double* b, double* x);

private:
    double dx; ///< Grid spacing in the x-direction.
    double dy; ///< Grid spacing in the y-direction.
    int Nx; ///< Total number of grid points in the x-direction.
    int Ny; ///< Total number of grid points in the y-direction.
    int Nx_local; ///< Number of local grid points in the x-direction.
    int Ny_local; ///< Number of local grid points in the y-direction.
    double* r; ///< Residual vector.
    double* k; ///< Direction vector.
    double* z; ///< Preconditioned residual vector.
    double* t; ///< Temporary vector for matrix-vector multiplication.

    const int p; ///< Number of processes in each grid direction.
    const int* coords; ///< Coordinates of the process in the grid.
    

    double* bufRecvTop; ///< Buffer for receiving data from the top neighbor.
    double* bufRecvBottom; ///< Buffer for receiving data from the bottom neighbor.
    double* bufRecvLeft; ///< Buffer for receiving data from the left neighbor.
    double* bufRecvRight; ///< Buffer for receiving data from the right neighbor.

    void ApplyOperator(double* k, double* t); ///< Applies the linear operator A to vector k, storing the result in t.
    void Precondition(double* k, double* t); ///< Applies the preconditioner to vector k, storing the result in t.
    void ImposeBC(double* k); ///< Imposes boundary conditions on vector k.
    void InitialiseBuffers(); ///< Initializes communication buffers for exchanging boundary data.
    void CleanUpBuffers(); ///< Cleans up dynamically allocated memory for communication buffers.
    
    Communicator* communicator; ///< Pointer to the communicator for parallel operations.

};
