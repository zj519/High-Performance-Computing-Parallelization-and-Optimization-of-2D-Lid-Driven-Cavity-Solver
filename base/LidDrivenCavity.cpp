#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
using namespace std;

#include <cblas.h>

///the macro IDX(I,J) calculates the 1D index of the array element based 
///on its row (J) and column (I) indices
#define IDX(I,J) ((J)*Nx + (I))

#include "LidDrivenCavity.h"
#include "SolverCG.h"

/**
 * @brief Constructor for LidDrivenCavity class.
 */
LidDrivenCavity::LidDrivenCavity()
{
}

/**
 * @brief Destructor for LidDrivenCavity class.
 */
LidDrivenCavity::~LidDrivenCavity()
{
    CleanUp();
}

/**
 * @brief Set the size of the computational domain.
 *
 * @param xlen Length of the domain in the x.
 * @param ylen Length of the domain in the y.
 */
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    this->Lx = xlen;
    this->Ly = ylen;
    UpdateDxDy();
}

/**
 * @brief Set the grid size for the computational domain.
 *
 * @param nx Number of grid points in the x.
 * @param ny Number of grid points in the y.
 */
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    this->Nx = nx;
    this->Ny = ny;
    UpdateDxDy();
}

/**
 * @brief Set the time step
 *
 * @param deltat The time step size.
 */
void LidDrivenCavity::SetTimeStep(double deltat)
{
    this->dt = deltat;
}

/**
 * @brief Set the final time
 *
 * @param finalt The end of simulation time.
 */
void LidDrivenCavity::SetFinalTime(double finalt)
{
    this->T = finalt;
}

/**
 * @brief Set the Reynolds number and calculate kinematic viscosity.
 *
 * @param re The Reynolds number to be set.
 */
void LidDrivenCavity::SetReynoldsNumber(double re)
{
    this->Re = re;
    this->nu = 1.0/re;  ///calculate kinematic viscosity
}

/**
 * @brief Initialize the LidDrivenCavity solver.
 */
void LidDrivenCavity::Initialise()
{
    CleanUp();

    v   = new double[Npts]();   ///initialize velocity array with a length of number of grid points
    s   = new double[Npts]();   ///initialize stream function array with a length of number of grid points
    tmp = new double[Npts]();   ///initialize temporary array with a length of number of grid points      
    cg  = new SolverCG(Nx, Ny, dx, dy); ///dynamically allocate memory for the new "SolverCG" object just created
}

/**
 * @brief Integrate the LidDrivenCavity solver over time.
 */
void LidDrivenCavity::Integrate()
{
    int NSteps = ceil(T/dt);    ///calculate number of timesteps
    for (int t = 0; t < NSteps; ++t)
    {   ///for loo to print out the step number and current time 
        std::cout << "Step: " << setw(8) << t
                  << "  Time: " << setw(8) << t*dt
                  << std::endl;
        Advance();  ///advance the simulation by one time step
    }
}

/**
 * @brief Write the solution to a file
 * @param file Name of the output file
 */
void LidDrivenCavity::WriteSolution(std::string file)
{
    double* u0 = new double[Nx*Ny]();   ///initialize array with size of Nx*Ny to store velocity component, u
    double* u1 = new double[Nx*Ny]();   ///initialize array with size of Nx*Ny to store velocity component, v
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u0[IDX(i,j)] =  (s[IDX(i,j+1)] - s[IDX(i,j)]) / dy; ///forward euler to calculate velocity component u from stream function
            u1[IDX(i,j)] = -(s[IDX(i+1,j)] - s[IDX(i,j)]) / dx; ///forward euler to calculate velocity component v from stream functio
        }
    }
    ///set boundary condition for the top such that u=U
    for (int i = 0; i < Nx; ++i) {
        u0[IDX(i,Ny-1)] = U;
    }

    std::ofstream f(file.c_str());
    std::cout << "Writing file " << file << std::endl;
    int k = 0;
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            k = IDX(i, j);  ///calculate index 'k' of the current grid point using IDX macro
            ///writes the following into a file: 
            ///x-cooridnate: i*dx; y-coordinate: j*dy
            ///v[k]: vorticity;    s[k]: stream function
            ///u0[k]:velocity component u;  u1[k]:velocity component v
            f << i * dx << " " << j * dy << " " << v[k] <<  " " << s[k] 
              << " " << u0[k] << " " << u1[k] << std::endl;
        }
        f << std::endl;
    }
    f.close();
    ///dealocate memory for arrays
    delete[] u0;
    delete[] u1;
}

/**
 * @brief Print the configuration of the LidDrivenCavity simulation
 * This function prints out config of the LidDrivenCavity simulation, including
 * grid size, spacing, domain length in x and y;
 * total number of grid points, size of time step, total number of time steps, 
 * input Reynolds number, 
 */
void LidDrivenCavity::PrintConfiguration()
{
    cout << "Grid size: " << Nx << " x " << Ny << endl;
    cout << "Spacing:   " << dx << " x " << dy << endl;
    cout << "Length:    " << Lx << " x " << Ly << endl;
    cout << "Grid pts:  " << Npts << endl;
    cout << "Timestep:  " << dt << endl;
    cout << "Steps:     " << ceil(T/dt) << endl;
    cout << "Reynolds number: " << Re << endl;
    cout << "Linear solver: preconditioned conjugate gradient" << endl;
    cout << endl;
    ///check if the time-step restriction is satisfies, it should <  0.25
    ///exit with -1 status if restriction not satisfied
    if (nu * dt / dx / dy > 0.25) {
        cout << "ERROR: Time-step restriction not satisfied!" << endl;
        cout << "Maximum time-step is " << 0.25 * dx * dy / nu << endl;
        exit(-1);
    }
}

/**
 * @brief clean up resources allocated by the simulation
 */
void LidDrivenCavity::CleanUp()
{
    if (v) {
        delete[] v;
        delete[] s;
        delete[] tmp;
        delete cg;
    }
}

/**
 * @brief calculate grid size in x and y, also total number of grid points
 */
void LidDrivenCavity::UpdateDxDy()
{
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    Npts = Nx * Ny;
}

/**
 * @brief the main algorithm to calculate vorticity and stream function based on 2.1 Algorithm
 */
void LidDrivenCavity::Advance()
{
    double dxi  = 1.0/dx;
    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;    ///calculate 1/(dx^2)
    double dy2i = 1.0/dy/dy;    ///calculate 1/(dy^2)

    // Boundary node vorticity
    for (int i = 1; i < Nx-1; ++i) {
        // bottom Equation 7
        v[IDX(i,0)]    = 2.0 * dy2i * (s[IDX(i,0)]    - s[IDX(i,1)]);
        // top  Equation 6
        v[IDX(i,Ny-1)] = 2.0 * dy2i * (s[IDX(i,Ny-1)] - s[IDX(i,Ny-2)])
                       - 2.0 * dyi*U;
    }
    for (int j = 1; j < Ny-1; ++j) {
        // left Equation 8
        v[IDX(0,j)]    = 2.0 * dx2i * (s[IDX(0,j)]    - s[IDX(1,j)]);
        // right Equation 9
        v[IDX(Nx-1,j)] = 2.0 * dx2i * (s[IDX(Nx-1,j)] - s[IDX(Nx-2,j)]);
    }

    // Compute interior vorticity Equation 10
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            v[IDX(i,j)] = dx2i*(
                    2.0 * s[IDX(i,j)] - s[IDX(i+1,j)] - s[IDX(i-1,j)])
                        + 1.0/dy/dy*(
                    2.0 * s[IDX(i,j)] - s[IDX(i,j+1)] - s[IDX(i,j-1)]);
        }
    }

    // Time advance vorticity   Equation 11
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            v[IDX(i,j)] = v[IDX(i,j)] + dt*(
                ( (s[IDX(i+1,j)] - s[IDX(i-1,j)]) * 0.5 * dxi
                 *(v[IDX(i,j+1)] - v[IDX(i,j-1)]) * 0.5 * dyi)
              - ( (s[IDX(i,j+1)] - s[IDX(i,j-1)]) * 0.5 * dyi
                 *(v[IDX(i+1,j)] - v[IDX(i-1,j)]) * 0.5 * dxi)
              + nu * (v[IDX(i+1,j)] - 2.0 * v[IDX(i,j)] + v[IDX(i-1,j)])*dx2i
              + nu * (v[IDX(i,j+1)] - 2.0 * v[IDX(i,j)] + v[IDX(i,j-1)])*dy2i);
        }
    }

    // Sinusoidal test case with analytical solution, which can be used to test
    // the Poisson solver
    /*
    const int k = 3;
    const int l = 3;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            v[IDX(i,j)] = -M_PI * M_PI * (k * k + l * l)
                                       * sin(M_PI * k * i * dx)
                                       * sin(M_PI * l * j * dy);
        }
    }
    */

    // Solve Poisson problem
    cg->Solve(v, s);
}