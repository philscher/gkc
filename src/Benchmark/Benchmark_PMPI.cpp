#include "mpi.h"
#include <iostream>

#include "Benchmark_PMPI.h"
#include <cstring>
static int nsend = 0;
  
struct Benchmark_PMPI::pmpi_vec gl_pmpi_vec;
 
std::vector<Benchmark_PMPI::pmpi_vec> Benchmark_PMPI::time_trace;
double Benchmark_PMPI::time_start;



Benchmark_PMPI::Benchmark_PMPI(Setup *setup, Parallel *_parallel, FileIO *_fileIO) : parallel(_parallel), fileIO(_fileIO)
{
   initData();

   time_start = MPI_Wtime();
}

 
Benchmark_PMPI::~Benchmark_PMPI()
{

  //// Write data 
  hid_t pmpiGroupTT = fileIO->newGroup("TimeTrace", pmpiGroup);
  
  //////////////////////// Set HDF-5 table for rank, Coordinates(6)
  hid_t str_tid = H5Tcopy(H5T_C_S1); H5Tset_size(str_tid, 16); H5Tset_strpad(str_tid, H5T_STR_NULLTERM);
  size_t rc_offset[]     = { HOFFSET(pmpi_vec, t), HOFFSET(pmpi_vec, dt), HOFFSET(pmpi_vec, dest), HOFFSET(pmpi_vec, bytes), 
                             HOFFSET(pmpi_vec, cmd) };
  size_t rc_sizes[]      = { sizeof(double), sizeof(double), sizeof(int), sizeof(int), sizeof(char) * 16 };
  hid_t rc_types[]       = { H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INT, H5T_NATIVE_INT, str_tid };
  const char *rc_names[] = { "Time", "dt", "Destination", "Bytes", "Command" }; 
 
for(int n = 0; n < parallel->numProcesses; n++) {
  //use auto-ptr
  TableAttr *eventTable = new TableAttr(pmpiGroupTT, "Rank_" + Setup::num2str(n), 5, rc_names, rc_offset, rc_types, rc_sizes, &_pmpi_vec);
  
  // Oh lord, this is so disgusting .... , check HDF-5 forum for non-collective HDF-5 table calls 
  int vecsize = parallel->reduce( n == parallel->myRank ? (int) Benchmark_PMPI::time_trace.size() : 0, Op::sum);
  char val[sizeof(pmpi_vec) * vecsize]; val[:] = 0.; 
  parallel->reduce( n == parallel->myRank ? (char *) Benchmark_PMPI::time_trace.data() : &val[0], Op::sum, DIR_ALL, (int) (sizeof(_pmpi_vec) * vecsize));
  eventTable->append((pmpi_vec *) &val[0], vecsize);

  delete eventTable;
}
  //if(n == parallel->myRank) eventTable->append(Benchmark_PMPI::time_trace.data(),  Benchmark_PMPI::time_trace.size());
  //TableAttr *eventTable = new TableAttr(pmpiGroupTT, "Rank_" + Setup::num2str(parallel->myRank), 2, rc_names, rc_offset, rc_types, rc_sizes, &_pmpi_vec);

  H5Gclose(pmpiGroupTT);
  H5Gclose(pmpiGroup);

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



void Benchmark_PMPI::initData()
{
  /////////// Setup Table with rank -vs coordinate
  pmpiGroup = fileIO->newGroup("/PMPI");
  
  //////////////////////// Set HDF-5 table for rank, Coordinates(6)
  struct rc { int mpi_rank; int coord[6]; } _rc;
  hsize_t   ull6 = 6ull; 
  hid_t coord_i6 = H5Tarray_create(H5T_NATIVE_INT, 1, &ull6);

  size_t rc_offset[]     = { HOFFSET(rc,mpi_rank), HOFFSET(rc, coord[0]) };
  size_t rc_sizes[]      = { sizeof(int), sizeof(int[6]) };
  hid_t rc_types[]       = { H5T_NATIVE_INT, coord_i6 };
  const char *rc_names[] = { "Rank", "Coordinate" }; 
   
  TableAttr *coordTable = new TableAttr(pmpiGroup, "Coordinates", 2, rc_names, rc_offset, rc_types, rc_sizes, &_rc);

  ////// Fill table
  int values[parallel->numProcesses][7]; values[:][:] = 0;

  values[parallel->myRank][0  ] = parallel->myRank;
  values[parallel->myRank][1:6] = parallel->Coord[0:6];

  parallel->reduce(&values[0][0], Op::sum, DIR_ALL, parallel->numProcesses * 7);
  coordTable->append(values, parallel->numProcesses);

  delete coordTable;
}

/*

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
  //return PMPI_Isend(start, count, datatype, dest, tag, comm, req);
  double call_time_start = MPI_Wtime();
  int ret =  PMPI_Isend(start, count, datatype, dest, tag, comm, req);
  double call_time_end   = MPI_Wtime();
 
  gl_pmpi_vec.t     = call_time_start - Benchmark_PMPI::time_start;
  gl_pmpi_vec.dt    = call_time_end   - call_time_start;
  gl_pmpi_vec.dest  = dest; 
  gl_pmpi_vec.bytes = count * sizeof(datatype); 
  std::strncpy(gl_pmpi_vec.cmd, "MPI_Isend", 9);
      
  Benchmark_PMPI::time_trace.push_back(gl_pmpi_vec);
  return ret;
}


} // extern "C" 
  */
