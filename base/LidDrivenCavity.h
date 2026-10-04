/**
 * @file LidDrivencavity.h
 * @brief Define the LidDrivenCavity class
 */

///prevent header file from being included multiple times in the same translation unit
#pragma once

#include <string>
using namespace std;

class SolverCG;

class LidDrivenCavity
{
public:
    /**
     * @brief Construct a new LidDrivenCavity object
     */
    LidDrivenCavity();
    
    /**
     * @brief Destroys the LidDrivenCavity object.
     */
    ~LidDrivenCavity();
    
    /**
     * @brief Sets the domain size of the cavity.
     * 
     * @param xlen The length of the domain in the x-direction.
     * @param ylen The length of the domain in the y-direction.
     */
    void SetDomainSize(double xlen, double ylen);
    
    /**
     * @brief Sets the grid size for discretizing the domain.
     * 
     * @param nx The number of grid points in the x-direction.
     * @param ny The number of grid points in the y-direction.
     */
    void SetGridSize(int nx, int ny);
    
    /**
     * @brief Sets the time step for time integration.
     * 
     * @param deltat The time step size.
     */
    void SetTimeStep(double deltat);
    
    /**
     * @brief Sets the final time for simulation.
     * 
     * @param finalt The final time for simulation.
     */
    void SetFinalTime(double finalt);
    
     /**
     * @brief Sets the Reynolds number for the flow.
     * 
     * @param Re The Reynolds number.
     */   
    void SetReynoldsNumber(double Re);

     /**
     * @brief Initializes the solver.
     */
    void Initialise();
    
    /**
     * @brief Integrates the flow equations over time.
     */
    void Integrate();
    
    /**
     * @brief Writes the solution data to a file.
     * 
     * Writes the solution data, including velocity components and pressure
     * @param file The name of the file to write the solution data to.
     */
    void WriteSolution(std::string file);
    
    /**
     * @brief Prints the configuration of the solver.
     */
    void PrintConfiguration();

public:
    double* v   = nullptr;  ///initialize velocity component array, v to a null pointer
    double* s   = nullptr;  ///initialize stream function array to a null pointer
    double* tmp = nullptr;  ///initialize a temporary array to a null pointer

    double dt   = 0.01; ///set time step size equals 0.01
    double T    = 1.0;  ///set final time equals 1.0
    double dx;  ///initialize grid spacing in x-direction
    double dy;  ///initialize grid spacing in y-direction
    int    Nx   = 9;  ///set number of grid points in x 
    int    Ny   = 9;  ///set number of grid poinys in y
    int    Npts = 81;  ///set total number of grid points
    double Lx   = 1.0;  ///set length of the domain in x
    double Ly   = 1.0;  ///set length of the domain in y
    double Re   = 10;  ///set Reynolds number 
    double U    = 1.0;  ///set the bulk velocity
    double nu   = 0.1;  ///set kinematic viscosity 

    SolverCG* cg = nullptr;  ///initialize conjugate gradient solver object to a null pointer
    
    ///clean up memory allocated for the solver
    void CleanUp();
    ///update the grid spacing
    void UpdateDxDy();
    ///advance the solution in time
    void Advance();
};