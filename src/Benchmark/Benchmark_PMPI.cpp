#include "mpi.h"

#include <iostream>

static int nsend = 0;
/*
extern "C" {
//int PMPI_send(void *start, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm );

int MPI_Send( void *start, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm )
{
 nsend++;
 std::cout << "Hi";
    PMPI_Send(start, count, datatype, dest, tag, comm);
 return 0;
}

int MPI_Isend(void*start, int count , MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *req) 
{
  std::cout << "MPI_ISend" << "d : " << dest << " t : " << tag << std::endl;

  return PMPI_Isend(start, count, datatype, dest, tag, comm, req);
}

}
*/
