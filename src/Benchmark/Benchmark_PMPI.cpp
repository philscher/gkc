#include "mpi.h"
#include <iostream>

#include "Benchmark_PMPI.h"

static int nsend = 0;


Benchmark_PMPI::Benchmark_PMPI(Setup *setup, Parallel *_parallel, FileIO *fileIO) : parallel(_parallel)
{


   /////////// Setup Table with rank -vs coordinate
   hid_t pmpiGroup = fileIO->newGroup("/PMPI");
  
   //////////////////////// Set Table for Events
   struct rp { int mpi_rank; int coord[6]; } _rp;
   hsize_t   ull6 = 6ull; 
   hid_t coord_i6 = H5Tarray_create(H5T_NATIVE_INT, 1, &ull6);

   size_t rp_offset[]     = { HOFFSET(rp,mpi_rank), HOFFSET(rp, coord[0]) };
   size_t rp_sizes[]      = { sizeof(int), sizeof(int[6]) };
   hid_t rp_types[]       = { H5T_NATIVE_INT, coord_i6 };
   const char *rp_names[] =  { "Rank", "Coord" }; 
   
   TableAttr *eventTable = new TableAttr(pmpiGroup, "Coords", 2, rp_names, rp_offset, rp_types, rp_sizes, &_rp);

   // Fill table
   int values[parallel->numProcesses][7]; values[:][:] = 0;

   values[parallel->myRank][0  ] = parallel->myRank;
   values[parallel->myRank][1:6] = parallel->Coord[0:6];

   parallel->reduce(&values[0][0], Op::sum, DIR_ALL, parallel->numProcesses * 7);
   eventTable->append(values, parallel->numProcesses);

   delete eventTable;
   H5Gclose(pmpiGroup);
   // save coordinates
   

}



/* 
void Benchmark_PMPI::start(std::string id)
{
      nsend = 0;
      
}



void Benchmark_PMPI::stop()
{


}

 * */

//int PMPI_send(void *start, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm );


////////////////////////////// Overloaded (P)MPI calls //////////////////////////

extern "C" { // required to avoid name mangling

int MPI_Send( void *start, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm )
{
 nsend++;
    int ret = PMPI_Send(start, count, datatype, dest, tag, comm);
 return ret;
}

int MPI_Isend(void*start, int count , MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *req) 
{
  return PMPI_Isend(start, count, datatype, dest, tag, comm, req);
  double ts = MPI_Wtime();
  int ret =  PMPI_Isend(start, count, datatype, dest, tag, comm, req);
  double dt = MPI_Wtime() - ts;
  
  std::cout << "MPI_ISend" << "destingatio(rank) : " << dest << " count : " << count << dt << std::endl;

}


}
