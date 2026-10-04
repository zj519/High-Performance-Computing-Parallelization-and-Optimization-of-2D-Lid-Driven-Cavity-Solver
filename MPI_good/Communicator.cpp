/**
 * @file Communicator.h
 * @brief Handles MPI communication for grid partitioning in parallel computing environments.
 */

#include "Communicator.h"
#include "mpi.h"
#include <iostream>
using namespace std;
#define IDX(I, J) ((J) * Nx_local + (I))

/**
 * @class Communicator
 * @brief Manages MPI communication for data exchange between grid partitions.
 *
 * This class is responsible for setting up and managing the grid communicator,
 * and for facilitating the exchange of boundary data between adjacent partitions
 * of a distributed grid in a parallel computing environment.
 */
Communicator::Communicator(MPI_Comm comm, int Nx_local, int Ny_local, int x_ini, int end_x, int y_ini, int end_y, int* coords, int p)
    : cart_comm(comm), Nx_local(Nx_local), Ny_local(Ny_local), x_ini(x_ini), end_x(end_x), y_ini(y_ini), end_y(end_y), coords(coords), p(p)
{
    // Initialize send buffers for grid boundary communication
    this->bufSendTop = new double[Nx_local];
    this->bufSendBottom = new double[Nx_local];
    this->bufSendLeft = new double[Ny_local];
    this->bufSendRight = new double[Ny_local];
    
}

/**
 * @brief Destructor to clean up allocated resources.
 */
Communicator::~Communicator()
{
    delete[] bufSendTop;
    delete[] bufSendBottom;
    delete[] bufSendLeft;
    delete[] bufSendRight;
}

/**
 * @brief Exchanges boundary data between adjacent partitions of the grid.
 *
 * This method performs non-blocking sends and receives to exchange boundary
 * data between the current partition and its adjacent top, bottom, left, and
 * right neighbors in the grid.
 *
 * @param Array The array containing the data to be sent to adjacent partitions.
 * @param receiveBufferTop Buffer to receive data from the top neighbor.
 * @param receiveBufferBottom Buffer to receive data from the bottom neighbor.
 * @param receiveBufferLeft Buffer to receive data from the left neighbor.
 * @param receiveBufferRight Buffer to receive data from the right neighbor.
 */
void Communicator::ExchangeBoundaryData(const double* Array, double* receiveBufferTop, double* receiveBufferBottom, double* receiveBufferLeft, double* receiveBufferRight) {
    MPI_Request requests[8]; // Send and receive in 4 directions, 2x4=8
    int requestCount = 0;

    int top, bottom, left, right;
    MPI_Cart_shift(cart_comm, 0, 1, &top, &bottom);
    MPI_Cart_shift(cart_comm, 1, 1, &left, &right);

    // Communicate with adjacent top neighbour if exists
    if (top != MPI_PROC_NULL) {
        MPI_Irecv(receiveBufferTop, Nx_local, MPI_DOUBLE, top, 0, cart_comm, &requests[requestCount++]);
        MPI_Isend(&Array[IDX(0, Ny_local - 1)], Nx_local, MPI_DOUBLE, top, 1, cart_comm, &requests[requestCount++]);
    }
    
    // communicate with adjacent bottom neighbour if exists
    if (bottom != MPI_PROC_NULL) {
        MPI_Irecv(receiveBufferBottom, Nx_local, MPI_DOUBLE, bottom, 1, cart_comm, &requests[requestCount++]);
        MPI_Isend(&Array[IDX(0, 0)], Nx_local, MPI_DOUBLE, bottom, 0, cart_comm, &requests[requestCount++]);
    }
    // communicate with adjacent left neighbour if exists
    if (left != MPI_PROC_NULL) {
        MPI_Irecv(receiveBufferLeft, Ny_local, MPI_DOUBLE, left, 2, cart_comm, &requests[requestCount++]);
        MPI_Isend(&Array[IDX(0, 0)], Ny_local, MPI_DOUBLE, left, 3, cart_comm, &requests[requestCount++]);
    }
    // communicate with adjacent right neighbour if exists
    if (right != MPI_PROC_NULL) {
        MPI_Irecv(receiveBufferRight, Ny_local, MPI_DOUBLE, right, 3, cart_comm, &requests[requestCount++]);
        MPI_Isend(&Array[IDX(Nx_local - 1, 0)], Ny_local, MPI_DOUBLE, right, 2, cart_comm, &requests[requestCount++]);
    }

    MPI_Waitall(requestCount, requests, MPI_STATUSES_IGNORE);
}

