#pragma once
/**
 * @class SolverCG
 * @brief A class for solving linear systems using the Conjugate Gradient method.
 */
class SolverCG
{
public:
     /**
     * @brief Constructor for the SolverCG class.
     * @param pNx Number of grid points in the x.
     * @param pNy Number of grid points in the y.
     * @param pdx Grid spacing in the x.
     * @param pdy Grid spacing in the y.
     * @param rank MPI rank of the current process.
     * @param numProcs Total number of MPI processes.
     */
    SolverCG(int pNx, int pNy, double pdx, double pdy);
    
    /**
     * @brief Destructor for the SolverCG class.
     */
    ~SolverCG();
    
    /**
     * @brief Solves the linear system Ax = b using the Conjugate Gradient method.
     * @param b Right-hand side vector.
     * @param x Solution vector.
     */
    void Solve(double* b, double* x);

public:
    double dx;
    double dy;
    int Nx;
    int Ny;
    double* r;
    double* p;
    double* z;
    double* t;
    int rank;
    int numProcs;
    int local_Nx;
    int local_Ny;

    void ApplyOperator(double* p, double* t);
    void Precondition(double* p, double* t);
    void ImposeBC(double* p);

};

