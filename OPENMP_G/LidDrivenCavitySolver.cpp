/**
 * @file main.cpp
 * @brief Entry point for the 2D lid-driven cavity solver using MPI for parallelization.
 * 
 * This program solves the incompressible Navier-Stokes equations in a square
 * domain using the finite difference method, parallelized with MPI. It accepts
 * various command-line options to specify the simulation parameters.
 */
#include <iostream>
#include "LidDrivenCavity.h"
#include "mpi.h"
#include <cmath> 
#include <omp.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/**
 * @brief Main function. Parses command-line arguments and runs the simulation.
 * 
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Returns 0 on successful completion, or 1 if an error occurs.
 */
int main(int argc, char **argv)
{    
   po::options_description opts(
        "Solver for the 2D lid-driven cavity incompressible flow problem");
    opts.add_options()
        ("Lx",  po::value<double>()->default_value(1.0),
                 "Length of the domain in the x-direction.")
        ("Ly",  po::value<double>()->default_value(1.0),
                 "Length of the domain in the y-direction.")
        ("Nx",  po::value<int>()->default_value(101),
                 "Number of grid points in x-direction.")
        ("Ny",  po::value<int>()->default_value(101),
                 "Number of grid points in y-direction.")
        ("dt",  po::value<double>()->default_value(0.005),
                 "Time step size.")
        ("T",   po::value<double>()->default_value(50),
                 "Final time.")
        ("Re",  po::value<double>()->default_value(1000),
                 "Reynolds number.")
        ("verbose",    "Be more verbose.")
        ("help",       "Print help message.");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << opts << endl;
        return 0;
    }
    
    /// Initialize MPI
    MPI_Init(&argc, &argv);
    // OpenMP thread initialization (optional, for demonstration)
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        if (threadID == 0) { // Master thread
            int nThreads = omp_get_num_threads();
            std::cout << "Running with " << nThreads << " OpenMP threads." << std::endl;
        }
    }
    
    double t0=MPI_Wtime();
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p=sqrt(size);
    if (p*p!=size) {
        if (rank == 0) {
            cout << "Errors: number of processors must be a perfect square" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    /// Create a Cartesian topology
    int dims[2] = {p, p};  
    int periods[2] = {0, 0};
    int coords[2];
    int reorder=0;
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(vm["Lx"].as<double>(), vm["Ly"].as<double>());
    solver->SetGridSize(vm["Nx"].as<int>(),vm["Ny"].as<int>());
    solver->SetTimeStep(vm["dt"].as<double>());
    solver->SetFinalTime(vm["T"].as<double>());
    solver->SetReynoldsNumber(vm["Re"].as<double>());
    solver->SetLocal(vm["Nx"].as<int>(), vm["Ny"].as<int>(), p, coords);
    solver->PrintConfiguration();
    solver->Initialise(cart_comm, coords, p);
    solver->WriteSolution("ic.txt");
    solver->Integrate();
    
    solver->WriteSolution("final.txt");
    double t_end=MPI_Wtime();
    double t_total=t_end-t0;
    cout << "Total execution time: " << t_total << " seconds." << endl;
    
    // clean up and Finalise MPI
    delete solver; 
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();

    return 0;
}