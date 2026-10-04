#pragma once

#include <string>
#include <tuple>
#include "mpi.h"

class SolverCG;
class Communicator;

class LidDrivenCavity {
public:
    // Constructor now takes a SimulationConfig object
    LidDrivenCavity();
    ~LidDrivenCavity();

    // Initialization function might no longer need separate setters for each parameter
    void Initialise(MPI_Comm comm, int* coords, int p);
    void Integrate();
    void WriteSolution(const std::string file);
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetLocal(int Nx, int Ny, int p, int *coords);
    void PrintConfiguration();
    // Tuple function for fetching neighbor values
    std::tuple<double, double, double, double> FetchNeighborValues(int i, int j, double* array, double* bufRecvLeft, double* bufRecvRight, double* bufRecvTop, double* bufRecvBottom);

public:
    double* v = nullptr;
    double* vv = nullptr;
    double* s = nullptr;
    double* tmp = nullptr;
    double dt   = 0.01;
    double T    = 1.0;
    double dx;
    double dy;
    int    Nx   = 9;
    int    Ny   = 9;
    int    Npts = 81;
    double Lx   = 1.0;
    double Ly   = 1.0;
    double Re   = 10;
    double U    = 1.0;
    double nu   = 0.1;

    double rank = 0; 
    double size = 0; 
    int* coords = nullptr;
    int p = 1;
    SolverCG* cg = nullptr;
    Communicator* communicator = nullptr;
    int Nx_local = Nx;
    int Ny_local = Ny;
    int Npts_local = Npts;


    double dxi;
    double dyi;
    double dx2i;
    double dy2i;
    
    int x_0;
    int x_end;
    int y_0;
    int y_end;

    double* bufRecvTop_V = nullptr;
    double* bufRecvBottom_V = nullptr;
    double* bufRecvLeft_V = nullptr;
    double* bufRecvRight_V = nullptr;

    double* bufRecvTop_S = nullptr;
    double* bufRecvBottom_S = nullptr;
    double* bufRecvLeft_S = nullptr;
    double* bufRecvRight_S = nullptr;

    void CleanUp();
    void UpdateDxDy();
    void Advance();
    void BoundaryV();
    void NodeVorticity();
    void InteriorV(int x_coord_0, int x_coord_end, int y_coord_0, int y_coord_end);
    void TimeAdvanceV(int x_coord_0, int x_coord_end, int y_coord_0, int y_coord_end);
};
