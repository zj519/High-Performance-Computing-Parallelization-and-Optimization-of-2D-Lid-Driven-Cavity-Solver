/**
 * @file Communicator.h
 * @brief Defines the Communicator class for managing MPI communications.
 */
#pragma once
#include "mpi.h"

/**
 * @class Communicator
 * @brief Manages MPI communications for exchanging boundary data between adjacent processes.
 *
 * This class is responsible for setting up a grid communicator for parallel processing,
 * initializing send and receive buffers, and performing non-blocking communication
 * to exchange boundary data with neighboring processes in a distributed grid.
 */
 
class Communicator {
public:
 /**
     * @brief Constructs a new Communicator object.
     *
     * Initializes MPI communication settings for the grid, including the communicator,
     * local grid sizes, start and end indices for x and y directions, coordinates of
     * the process within the grid, and the total number of processes in each direction.
     *
     * @param comm The MPI communicator for the grid.
     * @param Nx_local Number of grid points in the local domain in the x-direction.
     * @param Ny_local Number of grid points in the local domain in the y-direction.
     * @param x_ini Start index in the x-direction.
     * @param end_x End index in the x-direction.
     * @param y_ini Start index in the y-direction.
     * @param end_y End index in the y-direction.
     * @param coords Coordinates of the process in the grid.
     * @param p Number of processes in each grid direction (assumed to be square).
     */
    Communicator(MPI_Comm comm, int Nx_local, int Ny_local, int x_ini, int end_x, int y_ini, int end_y, int* coords, int p);
    
    /**
     * @brief Destroys the Communicator object.
     *
     * Cleans up dynamically allocated memory for send buffers.
     */
     
    ~Communicator();
        /**
     * @brief Performs non-blocking communication to exchange boundary data with neighboring processes.
     *
     * Initiates non-blocking sends and receives to exchange boundary data of the specified
     * array with the top, bottom, left, and right neighbors, if they exist.
     *
     * @param Array The array containing the data to be sent to adjacent processes.
     * @param bufRecvTop Buffer for receiving data from the top neighbor.
     * @param bufRecvBottom Buffer for receiving data from the bottom neighbor.
     * @param bufRecvLeft Buffer for receiving data from the left neighbor.
     * @param bufRecvRight Buffer for receiving data from the right neighbor.
     */
    void ExchangeBoundaryData(const double* Array, double* bufRecvTop, double* bufRecvBottom, double* bufRecvLeft, double* bufRecvRight);

    const MPI_Comm cart_comm; ///< MPI communicator for the grid.
    const int Nx_local, Ny_local; ///< Local grid sizes in x and y directions.
    const int x_ini, end_x; ///< Start and end indices in the x-direction.
    const int y_ini, end_y; ///< Start and end indices in the y-direction.
    const int* coords; ///< Coordinates of the process in the grid.
    const int p; ///< Number of processes in each grid direction.

private:
    double* bufSendTop; ///< Buffer for sending data to the top neighbor.
    double* bufSendBottom; ///< Buffer for sending data to the bottom neighbor.
    double* bufSendLeft; ///< Buffer for sending data to the left neighbor.
    double* bufSendRight; ///< Buffer for sending data to the right neighbor.
};


