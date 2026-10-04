/**
 * @file LidDrivenCavity.h
 * @brief Header file for the LidDrivenCavity class.
 *
 * This file declares the LidDrivenCavity class and its member functions used
 * to solve the 2D lid-driven cavity problem using a finite difference method,
 * parallelized with MPI. The solver utilizes a Conjugate Gradient method for
 * solving the linear system of equations.
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <tuple>
#include <cblas.h>
#include <vector>
#define IDX(I, J) ((J) * Nx_local + (I))
#define IDX_glob(I, J) ((J) * Nx + (I))

#include "LidDrivenCavity.h"
#include "SolverCG.h"
#include "Communicator.h"
#include "mpi.h"
using namespace std;

/**
 * @class LidDrivenCavity
 * @brief Solves the 2D lid-driven cavity flow problem.
 *
 * The LidDrivenCavity class encapsulates the data and methods necessary to
 * simulate and solve the lid-driven cavity problem on a uniform grid. It
 * includes functionalities for setting up the simulation domain, initializing
 * the flow field, applying boundary conditions, advancing the solution in time,
 * and writing the results to a file.
 */
LidDrivenCavity::LidDrivenCavity()
{
   
}

LidDrivenCavity::~LidDrivenCavity()
{
    CleanUp();
}


/**
 * @fn void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
 * @brief Sets the physical dimensions of the simulation domain.
 * @param xlen Length of the domain in the x-direction.
 * @param ylen Length of the domain in the y-direction.
 */
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    this->Lx = xlen;
    this->Ly = ylen;
    UpdateDxDy();
}

/**
 * @fn void LidDrivenCavity::SetGridSize(int nx, int ny)
 * @brief Sets the size of the computational grid.
 * @param nx Number of grid points in the x-direction.
 * @param ny Number of grid points in the y-direction.
 */
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    this->Nx = nx;
    this->Ny = ny;
}

/**
 * @fn void LidDrivenCavity::SetTimeStep(double deltat)
 * @brief Sets the time step for the simulation.
 * @param deltat The time step size.
 */
void LidDrivenCavity::SetTimeStep(double deltat)
{
    this->dt = deltat;
}

/**
 * @fn void LidDrivenCavity::SetFinalTime(double finalt)
 * @brief Sets the final time for the simulation.
 * @param finalt The final simulation time.
 */
void LidDrivenCavity::SetFinalTime(double finalt)
{
    this->T = finalt;
}

/**
 * @fn void LidDrivenCavity::SetReynoldsNumber(double re)
 * @brief Sets the Reynolds number for the simulation.
 * @param re The Reynolds number.
 */

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    this->Re = re;
    this->nu = 1.0 / re;
}

/**
 * @fn void LidDrivenCavity::SetLocal(int Nx, int Ny, int p, int* coords)
 * @brief Sets local variables and computes the start, end indices for the local sub grid.
 * @param Nx Total number of grid points in x-direction.
 * @param Ny Total number of grid points in y-direction.
 * @param p The number of processes in one dimension.
 * @param coords The coordinates of the current process in the process grid.
 */
void LidDrivenCavity::SetLocal(int Nx, int Ny, int p, int* coords) {
    // Lambda function to calculate start, end indices and local grid size
    auto calculate_indices = [p](int N, int coord, int& start, int& end, int& N_local) {
        int extra_points = N % p;
        int min_points = N / p;
        int remainder = (coord < extra_points) ? 1 : 0;
        
        start = coord * (min_points + remainder);
        if (remainder == 1) {
            // Process gets an extra point
            end = start + min_points + 1;
        } else {
            // Process gets the minimum number of points, adjust start for processes after the extra points
            start += extra_points;
            end = start + min_points;
        }
        N_local = end - start;
    };
    // Calculate indices for X direction
    calculate_indices(Nx, coords[1], this->x_0, this->x_end, this->Nx_local);
    // Calculate indices for Y direction
    calculate_indices(Ny, coords[0], this->y_0, this->y_end, this->Ny_local);

    UpdateDxDy(); 
}

/**
 * @fn void LidDrivenCavity::Initialise(MPI_Comm comm, int *coords, int p)
 * @brief Initializes the simulation, setting up MPI communicators and allocating buffers.
 * @param comm The MPI communicator for the Cartesian topology.
 * @param coords The coordinates of the current process in the Cartesian grid.
 * @param p The number of processes per dimension in the Cartesian grid.
 */
void LidDrivenCavity::Initialise(MPI_Comm comm, int *coords, int p)
{
    bufRecvTop_V = new double[Nx_local];
    bufRecvBottom_V = new double[Nx_local];
    bufRecvLeft_V = new double[Ny_local];
    bufRecvRight_V = new double[Ny_local];

    bufRecvTop_S = new double[Nx_local];
    bufRecvBottom_S = new double[Nx_local];
    bufRecvLeft_S = new double[Ny_local];
    bufRecvRight_S = new double[Ny_local];

    int x_ini = coords[1] == 0 ? 1 : 0;
    int end_x = coords[1] == p - 1 ? Nx_local - 1 : Nx_local;
    int y_ini = coords[0] == p - 1 ? 1 : 0;
    int end_y = coords[0] == 0 ? Ny_local - 1 : Ny_local;
    

    communicator = new Communicator(comm, Nx_local, Ny_local, x_ini, end_x, y_ini, end_y, coords, p);
    cg = new SolverCG(Nx_local, Ny_local, dx, dy, communicator);
    this->coords = coords;
    this->p = p;
    v = new double[Npts_local](); 
    vv = new double[Npts_local](); 
    s = new double[Npts_local](); 
    tmp = new double[Npts_local]();
    
}

/**
 * @fn void LidDrivenCavity::Integrate()
 * @brief Advances the solution in time until the final time is reached.
 */
void LidDrivenCavity::Integrate()
{
    int NSteps = ceil(T / dt);
    for (int t = 0; t < NSteps; ++t)
    {
        // std::cout << "Step: " << setw(8) << t
        //           << "  Time: " << setw(8) << t * dt
        //           << std::endl;
        Advance();
    }
}

/**
 * @fn void LidDrivenCavity::PrintConfiguration()
 * @brief Prints the current configuration of the simulation to standard output.
 */
void LidDrivenCavity::PrintConfiguration()
{
    cout << "Grid size: " << Nx << " x " << Ny << endl;
    cout << "Spacing:   " << dx << " x " << dy << endl;
    cout << "Length:    " << Lx << " x " << Ly << endl;
    cout << "Grid pts:  " << Npts << endl;
    cout << "Timestep:  " << dt << endl;
    cout << "Steps:     " << ceil(T / dt) << endl;
    cout << "Reynolds number: " << Re << endl;
    cout << "Linear solver: preconditioned conjugate gradient" << endl;
    cout << endl;
    if (nu * dt / dx / dy > 0.25)
    {
        cout << "ERROR: Time-step restriction not satisfied!" << endl;
        cout << "Maximum time-step is " << 0.25 * dx * dy * Re << endl; 
        exit(-1);
    }
}

/**
 * @fn void LidDrivenCavity::CleanUp()
 * @brief Cleans up resources used by the LidDrivenCavity instance.
 */
void LidDrivenCavity::CleanUp()
{
    if (v)
    {
        delete[] v;
        delete[] s;
        delete[] tmp;
        delete cg;
        delete communicator;
        delete[] vv;

    }
}

/**
 * @fn void LidDrivenCavity::UpdateDxDy()
 * @brief Updates the grid spacing based on the current domain size and grid size.
 */
void LidDrivenCavity::UpdateDxDy()
{
    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    Npts_local = Nx_local * Ny_local;
    Npts = Nx * Ny;
}

/**
 * @fn void LidDrivenCavity::Advance()
 * @brief Advances the flow field by one time step using the vorticity-streamfunction formulation.
 */
void LidDrivenCavity::Advance()
{
    this->dxi = 1.0 / dx;
    this->dyi = 1.0 / dy;
    this->dx2i = 1.0 / dx / dx;
    this->dy2i = 1.0 / dy / dy;
    int x_ini = communicator->x_ini;
    int end_x = communicator->end_x;
    int y_ini = communicator->y_ini;
    int end_y = communicator->end_y;

    BoundaryV();   
    InteriorV(x_ini, end_x, y_ini, end_y); 
    TimeAdvanceV(x_ini, end_x, y_ini, end_y); 
    cg->Solve(v, s);
}

/**
 * @fn void LidDrivenCavity::BoundaryV()
 * @brief Applies boundary conditions for the vorticity at the cavity walls.
 */
void LidDrivenCavity::BoundaryV() {
    // Boundary condition for vorticity at the walls, depends on the position of the process.
    if (coords[0] == 0 || coords[0] == p - 1) {
         #pragma omp parallel for
        for (int i = 0; i < Nx_local; ++i) {
            // Top boundary
            if (coords[0] == 0) {
                vv[IDX(i, Ny_local - 1)] = 2.0 * dy2i * (s[IDX(i, Ny_local - 1)] - s[IDX(i, Ny_local - 2)]) - 2.0 * dyi * U;
            }
            // Bottom boundary
            if (coords[0] == p - 1) {
                vv[IDX(i, 0)] = 2.0 * dy2i * (s[IDX(i, 0)] - s[IDX(i, 1)]);
            }
        }
    }
    if (coords[1] == 0 || coords[1] == p - 1) {
         #pragma omp parallel for
        for (int j = 0; j < Ny_local; ++j) {
            // Left boundary
            if (coords[1] == 0) {
                vv[IDX(0, j)] = 2.0 * dx2i * (s[IDX(0, j)] - s[IDX(1, j)]);
            }
            // Right boundary
            if (coords[1] == p - 1) {
                vv[IDX(Nx_local - 1, j)] = 2.0 * dx2i * (s[IDX(Nx_local - 1, j)] - s[IDX(Nx_local - 2, j)]);
            }
        }
    }
}

/**
 * @fn void LidDrivenCavity::InteriorV(int x_coord_0, int x_coord_end, int y_coord_0, int y_coord_end)
 * @brief Updates the vorticity in the interior of the domain.
 * @param x_coord_0 Starting x-coordinate for the update region.
 * @param x_coord_end Ending x-coordinate for the update region.
 * @param y_coord_0 Starting y-coordinate for the update region.
 * @param y_coord_end Ending y-coordinate for the update region.
 */
void LidDrivenCavity::InteriorV(int x_coord_0, int x_coord_end, int y_coord_0, int y_coord_end) {
    communicator->ExchangeBoundaryData(s, bufRecvTop_S, bufRecvBottom_S, bufRecvLeft_S, bufRecvRight_S);
     #pragma omp parallel for
    for (int j = y_coord_0; j < y_coord_end; ++j) {
        for (int i = x_coord_0; i < x_coord_end; ++i) {
            auto [adjacentLeft, adjacentRight, adjacentBottom, adjacentTop] = FetchNeighborValues(i, j, s, 
                        bufRecvLeft_S, bufRecvRight_S, bufRecvTop_S, bufRecvBottom_S);

        
            vv[IDX(i, j)] = dx2i * (2.0 * s[IDX(i, j)] - adjacentRight - adjacentLeft)
                            + dy2i * (2.0 * s[IDX(i, j)] - adjacentTop - adjacentBottom);
        }
    }
}


/**
 * @brief Advances the vorticity field in time, applying boundary conditions and updating interior points.
 * 
 * This function advances the vorticity field in time using the current stream function values. It applies
 * boundary conditions for vorticity based on the no-slip condition at walls and updates the interior vorticity
 * values by solving the vorticity transport equation. MPI communication is used to exchange boundary data with
 * neighboring processes.
 * 
 * @param x_coord_0 Starting x-index for processing the local domain.
 * @param x_coord_end Ending x-index for processing the local domain.
 * @param y_coord_0 Starting y-index for processing the local domain.
 * @param y_coord_end Ending y-index for processing the local domain.
 */
void LidDrivenCavity::TimeAdvanceV(int x_coord_0, int x_coord_end, int y_coord_0, int y_coord_end) {
    // Ensuring proper communication for boundary values
    communicator->ExchangeBoundaryData(vv, bufRecvTop_V, bufRecvBottom_V, bufRecvLeft_V, bufRecvRight_V);
    communicator->ExchangeBoundaryData(s, bufRecvTop_S, bufRecvBottom_S, bufRecvLeft_S, bufRecvRight_S);

    for (int j = y_coord_0; j < y_coord_end; ++j) {
        for (int i = x_coord_0; i < x_coord_end; ++i) {
            // Fetching neighbor values for stream function 's' and vorticity 'vv'
            auto [S_left, S_right, S_bottom, S_top] = FetchNeighborValues(i, j, s, bufRecvLeft_S, bufRecvRight_S, bufRecvTop_S, bufRecvBottom_S);
            auto [V_left, V_right, V_bottom, V_top] = FetchNeighborValues(i, j, vv, bufRecvLeft_V, bufRecvRight_V, bufRecvTop_V, bufRecvBottom_V);
            
            v[IDX(i, j)] =   vv[IDX(i,j)] + dt*(( (S_right- S_left) * 0.5 * dt*(V_top - V_bottom) * 0.5 * dyi)- ( (S_top - S_bottom) * 0.5 * dyi
                        *(V_right - V_left) * 0.5 * dxi)+ nu * (V_right - 2.0 * vv[IDX(i,j)] + V_left)*dx2i+ nu * (V_bottom - 2.0 *vv[IDX(i,j)] + V_top)*dy2i);
                                    
        }
    }
}

/**
 * @brief Fetches the neighboring values required for numerical schemes from the local and received data.
 * 
 * This function fetches the values of neighboring points required for discretization schemes, taking into account
 * the position of the process in the MPI process grid to determine if neighbors are within the local domain or in
 * the halo region received from adjacent processes.
 * 
 * @param i x-index of the point whose neighbors are to be fetched.
 * @param j y-index of the point whose neighbors are to be fetched.
 * @param array Pointer to the local array containing the current process's data.
 * @param bufRecvLeft Pointer to the received left boundary data.
 * @param bufRecvRight Pointer to the received right boundary data.
 * @param bufRecvTop Pointer to the received top boundary data.
 * @param bufRecvBottom Pointer to the received bottom boundary data.
 * @return A tuple containing the values of the left, right, bottom, and top neighbors.
 */
std::tuple<double, double, double, double> LidDrivenCavity::FetchNeighborValues(int i, int j, double* array, double* bufRecvLeft, double* bufRecvRight, double* bufRecvTop, double* bufRecvBottom) 
{
    double adjacentLeft = (i == 0) ? bufRecvLeft[j] : array[IDX(i - 1, j)];
    double adjacentRight = (i == Nx_local-1) ? bufRecvRight[j] : array[IDX(i + 1, j)];
    double adjacentBottom = (j == 0) ? bufRecvBottom[i] : array[IDX(i, j - 1)];
    double adjacentTop = (j == Ny_local-1) ? bufRecvTop[i] : array[IDX(i, j + 1)];

    return std::make_tuple(adjacentLeft, adjacentRight, adjacentBottom, adjacentTop);
}

/**
 * @brief Writes the simulation results to a file in a structured format.
 * 
 * This function writes the simulation results, including velocity components and vorticity at each grid point,
 * to a specified file. It first gathers all data to the root process using MPI reductions and then the root
 * process writes the data to the file. The function handles the conversion from local to global indices and
 * ensures that only the root process performs the file I/O operations.
 * 
 * @param file The name of the file where the results will be written.
 */
void LidDrivenCavity::WriteSolution(std::string file)
{   
    int rank;
    MPI_Comm_rank(communicator->cart_comm, &rank);
    communicator->ExchangeBoundaryData(s, bufRecvTop_S, bufRecvBottom_S, bufRecvLeft_S, bufRecvRight_S); 

    int x_ini = communicator->x_ini;
    int end_x =communicator->end_x;
    int y_ini = communicator->y_ini;
    int end_y = communicator->end_y;
    int Npts_local=Nx_local*Ny_local;
    double *u0 = new double[Npts_local]();
    double *u1 = new double[Npts_local]();

    for (int j = y_ini; j < end_y; ++j)
    {
        for (int i = x_ini; i < end_x; ++i)
        {   
            double adjacentRight = (coords[1] < p - 1 && i == Nx_local-1) ? bufRecvRight_S[j] : s[IDX(i + 1, j)];
            double adjacentTop = (coords[0] > 0 && j == Ny_local-1) ? bufRecvTop_S[i] : s[IDX(i, j + 1)];

            u0[IDX(i, j)] = (adjacentTop - s[IDX(i, j)]) / dy;
            u1[IDX(i, j)] = -(adjacentRight - s[IDX(i, j)]) / dx;
        }
    }

    if (coords[0] == 0){
        for (int i = 0; i < Nx_local; ++i)
        {
            u0[IDX(i, Ny_local - 1)] = U;
        }
    } 

    double* u0_out = new double[Npts]();
    double* u1_out = new double[Npts]();
    double* v_out = new double[Npts]();
    double* s_out = new double[Npts]();

    double* u0_glob = new double[Npts]();
    double* u1_glob = new double[Npts]();
    double* v_glob = new double[Npts]();
    double* s_glob = new double[Npts]();


    for (int j = y_0; j < y_end; ++j)
    {
        for (int i = x_0; i < x_end; ++i)
        {   
            v_out[IDX_glob(i, j)] = v[IDX(i - x_0, j - y_0)];
            s_out[IDX_glob(i, j)] = s[IDX(i - x_0, j - y_0)];
            u0_out[IDX_glob(i, j)] = u0[IDX(i - x_0, j - y_0)];
            u1_out[IDX_glob(i, j)] = u1[IDX(i - x_0, j - y_0)];
            
        }
    }
    
    //sum each individual to form the global output 

    MPI_Allreduce(u0_out, u0_glob, Npts, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
    MPI_Allreduce(u1_out, u1_glob, Npts, MPI_DOUBLE, MPI_SUM,communicator->cart_comm);
    MPI_Allreduce(v_out, v_glob, Npts, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
    MPI_Allreduce(s_out, s_glob, Npts, MPI_DOUBLE, MPI_SUM, communicator->cart_comm);
    if (rank == 0){
        std::ofstream f(file.c_str());
        std::cout << "Writing file " << file << std::endl;
        int k = 0;
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                k = IDX_glob(i, j);
                f << i * dx << " " << j * dy << " " << v_glob[k] << " " << s_glob[k]
                        << " " << u0_glob[k] << " " << u1_glob[k] << std::endl;

            }
            f << std::endl;
        }
        f.close();

        delete[] u0;
        delete[] u1;
        delete[] u0_out;
        delete[] u1_out;
        delete[] v_out;
        delete[] s_out;
        delete[] u0_glob;
        delete[] u1_glob;
        delete[] v_glob;
        delete[] s_glob;
    }

}