/*
 * =====================================================================================
 *
 *       Filename: Parallel.cpp
 *
 *    Description: Implementation of Parallelization options. Now mainly uses
 *                 OpenMPI
 *
 *         Author: Paul P. Hilscher (2010), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#include "Global.h"
#include "Parallel.h"
#include "FileIO.h"
#include "Tools/System.h"

#include <array>

// global process rank needed for output
int process_rank;
int threadID;

// MPI_Error handler so we can use check and backtrace
void check_mpi(MPI_Comm *comm, int *err_code, ...) 
{

  char string[MPI_MAX_ERROR_STRING]; 
  int len;
  
  MPI_Error_string(*err_code, string, &len);

  std::cerr << "MPI Error : " << string << std::endl;  
  System::printStackTrace();
  
  if (*err_code != MPI_SUCCESS) check(-1, DMESG(std::string(string)), true);

}


Parallel::Parallel(Setup *setup) 
{

  // initialize some basic parameter
  myRank = 0;  threadID = 0; numThreads = 1; numProcesses = 1, master_rank = 0; 
  useOpenMP = false; useMPI = true;
  Coord[:] = 0;

  decomposition[:] = 1;

#ifdef GKC_PARALLEL_OPENMP
  useOpenMP = true;
  #pragma omp parallel
  if(omp_get_thread_num() == 0) numThreads = omp_get_num_threads();

#endif // GKC_PARALLEL_OPENMP

  for(int d = DIR_X; d < DIR_SIZE; d++) Comm[d] = MPI_COMM_NULL;

  //////////////////////// Set Message tags, enumerate through ////////////////////////
  
  int i = 14665; // a unique (random) number 
  
  for(int dir = DIR_X; dir <= DIR_S; dir++) {

    Talk[dir].psf_msg_tag[0] = ++i; Talk[dir].psf_msg_tag[1] = ++i;
    Talk[dir].phi_msg_tag[0] = ++i; Talk[dir].phi_msg_tag[1] = ++i;
  }
   
  // Initialize MPI
  int provided = 0; 
  int required = numThreads == 1 ? MPI_THREAD_SINGLE : MPI_THREAD_SERIALIZED;
  
  MPI_Init_thread(&setup->argc, &setup->argv, required, &provided);
  //check(provided < required ? -1 : 1, DMESG("MPI : Thread level support not available (use only one thread)"));

  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  //////////////////////// Set decomposition //////////////////////////////////
  
  std::vector<std::string> decomp = Setup::split(setup->get("Parallel.Decomposition","Auto"), ":");

  if(decomp.size() > 6) check(-1, DMESG("Decomposition only up to six dimensions"));

  if(decomp[0] == "Auto")  getAutoDecomposition(numProcesses);
  else for(int dir = DIR_X; dir < decomp.size() && (dir <= DIR_S); dir++) decomposition[dir] = std::stoi(decomp[dir]);
 
  checkValidDecomposition(setup);
   
  /////////////////// Create 6-D Cartesian Grid ////////////////////
  int periods[6] = { 1, 1, 1, 0, 0, 0};
  // Remove OpenMP decomposition, for the MPI call 
  decomposition[DIR_Y] = 1; 
  MPI_Cart_create(MPI_COMM_WORLD, 6, decomposition, periods, 1, &Comm[DIR_ALL]);

  // Get Cart coordinates to set 
  MPI_Comm_rank(Comm[DIR_ALL],&myRank); 
  process_rank = myRank; // process rank is global, needed ?

  MPI_Cart_coords(Comm[DIR_ALL], myRank, 6, Coord);
 
  // Get combined positions 
  Coord[DIR_XYZ]   = ((Coord[DIR_X] == 0) && (Coord[DIR_Y] == 0) && (Coord[DIR_Z] == 0)) ? 0 : -1;
  Coord[DIR_VMS]   = ((Coord[DIR_V] == 0) && (Coord[DIR_M] == 0) && (Coord[DIR_S] == 0)) ? 0 : -1;
  Coord[DIR_VM ]   = ((Coord[DIR_V] == 0) && (Coord[DIR_M] == 0))                        ? 0 : -1;
  Coord[DIR_MS ]   = ((Coord[DIR_M] == 0) && (Coord[DIR_S] == 0))                        ? 0 : -1;

  ///////////////////  Set SubCommunicators for various directions //////////////////////////
  // Note , we do not need to translate ranks, because dirMaster is always used with appropriate Comm
  int coord_master[6] = { 0, 0, 0, 0, 0, 0 };
    
  //////////////////////////////////////////////////////////////////////////////
  //
  //  We set the subcommunictor in direction dir using remain_dims
  //
  //  We prefer MPI_Cart_sub over MPI_Comm_split as the former includes
  //  more information about the topology and may results in more efficient
  //  communicators. 
  //  better use intializer list
  //
  
  auto setSubCommunicator = [=](int dir, int remain_dims[6]) 
  {
    int coord_master[6] = { 0, 0, 0, 0, 0, 0 };
    MPI_Cart_sub (Comm[DIR_ALL], remain_dims, &Comm[dir]);
    MPI_Cart_rank(Comm[dir    ], coord_master, &dirMaster[dir]);
  };
  /*
    setSubCommunicator(DIR_X  , { true , false, false, false, false, false } );
    setSubCommunicator(DIR_Z  , { false, false, true , false, false, false } );         
    setSubCommunicator(DIR_V  , { false, false, false, true , false, false } );
    setSubCommunicator(DIR_M  , { false, false, false, false, true , false } );
    setSubCommunicator(DIR_S  , { false, false, false, false, false, true  } );
    setSubCommunicator(DIR_MS , { false, false, false, false, true , true  } );
    setSubCommunicator(DIR_VM , { false, false, false, true , true , false } );
    setSubCommunicator(DIR_VMS, { false, false, false, true , true , true  } );
  */
  
  // Communicator for X
  int remain_dim_X  [6] = { true , false, false, false, false, false }; setSubCommunicator(DIR_X  , remain_dim_X  );      
  int remain_dim_Z  [6] = { false, false, true , false, false, false }; setSubCommunicator(DIR_Z  , remain_dim_Z  );        
  int remain_dim_V  [6] = { false, false, false, true , false, false }; setSubCommunicator(DIR_V  , remain_dim_V  );
  int remain_dim_M  [6] = { false, false, false, false, true , false }; setSubCommunicator(DIR_M  , remain_dim_M  );
  int remain_dim_S  [6] = { false, false, false, false, false, true  }; setSubCommunicator(DIR_S  , remain_dim_S  );
  int remain_dim_MS [6] = { false, false, false, false, true , true  }; setSubCommunicator(DIR_MS , remain_dim_MS );
  int remain_dim_XYZ[6] = { true , true , true , false, false, false }; setSubCommunicator(DIR_XYZ, remain_dim_XYZ);
  int remain_dim_VMS[6] = { false, false, false, true , true , true  }; setSubCommunicator(DIR_VMS, remain_dim_VMS);
  int remain_dim_VM [6] = { false, false, false, true , true , false }; setSubCommunicator(DIR_VM , remain_dim_VM );
  int remain_dim_XYZVM[6] = { true, true, true , true , true , false }; setSubCommunicator(DIR_XYZVM, remain_dim_XYZVM);
    
  ///////////////////////////////////////////////////////////////////////
  //
  //   Get ranks of neigbouring processes for Vlasov equation boundary 
  //   exchange, e.g. in X
  //
  //   ...|(x=n-2)|(x=n-1) rank_l|(x=n) my_rank|(x=n+1) rank_u|(x=+2)| ...
  //
  //   Note that, x,z are periodic, thus (n=-1) corresponds to (n=Nx-1)
  //
  auto setNeighbourRank = [=](int dir) 
  {
    int rank = myRank; // is this necessary ?!
    MPI_Cart_shift(Comm[DIR_ALL], dir, 1, &rank, &Talk[dir].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], dir,-1, &rank, &Talk[dir].rank_l);
  };
  //for(int dir: { DIR_X, DIR_Z, DIR_V, DIR_M, DIR_S }) setNeighbourRank(dir);
  int dirs[] = { DIR_X, DIR_Z, DIR_V, DIR_M, DIR_S };
  for(int i = 0; i < 5; i++) setNeighbourRank(dirs[i]);
  
  /////////////////// Setup own error handle,  this needed to enable backtracing when debugging
  MPI_Errhandler my_errhandler;
  MPI_Errhandler_create(&check_mpi, &my_errhandler);
  MPI_Errhandler_set(Comm[DIR_ALL], my_errhandler);

}

Parallel::~Parallel() 
{
  // free MPI communicators (Comm[] was initialized with MPI_COMM_NULL)
  for(int n=DIR_X; n<DIR_SIZE; n++) if(Comm[n] != MPI_COMM_NULL) MPI_Comm_free(&Comm[n]);

  check(MPI_Finalize(), DMESG("MPI Error"));
 
  return;
}


void Parallel::updateBoundaryVlasov(CComplex *Sendu, CComplex *Sendl, CComplex *Recvu, CComplex  *Recvl, int num, int dir)
{
  
 // OPTIM : If we don't decompose,  copy Sendu->Recvl, .., without going through MPI

#ifdef GKC_PARALLEL_MPI
     
  auto mpi_type = getMPIDataType(typeid(CComplex)); 

  bool nonblocking = true;
   
  if(nonblocking) {
     
    MPI_Irecv(Recvl, num, mpi_type, Talk[dir].rank_l, Talk[dir].psf_msg_tag[0], Comm[DIR_ALL], &Talk[dir].psf_msg_req[0]);
    MPI_Isend(Sendu, num, mpi_type, Talk[dir].rank_u, Talk[dir].psf_msg_tag[0], Comm[DIR_ALL], &Talk[dir].psf_msg_req[1]);
    MPI_Irecv(Recvu, num, mpi_type, Talk[dir].rank_u, Talk[dir].psf_msg_tag[1], Comm[DIR_ALL], &Talk[dir].psf_msg_req[2]); 
    MPI_Isend(Sendl, num, mpi_type, Talk[dir].rank_l, Talk[dir].psf_msg_tag[1], Comm[DIR_ALL], &Talk[dir].psf_msg_req[3]);
      
    if(Talk[dir].rank_u == MPI_PROC_NULL)   Recvl[0:num] = 0.e0;
    if(Talk[dir].rank_l == MPI_PROC_NULL)   Recvu[0:num] = 0.e0;

  } else {

    MPI_Status status[2];
    if(Talk[dir].rank_u == MPI_PROC_NULL)   std::cout << "I am Null ! rank_u" << std::endl << std::flush;
    if(Talk[dir].rank_l == MPI_PROC_NULL)   std::cout << "I am Null ! rank_l" << std::endl << std::flush;

    int tag[2] = { 123123, 123135 }; 
    MPI_Sendrecv(Recvl, num, mpi_type, Talk[dir].rank_u, tag[0], Recvl, num, mpi_type, Talk[dir].rank_l, tag[0], Comm[DIR_ALL], &status[0]);
    MPI_Sendrecv(Recvu, num, mpi_type, Talk[dir].rank_l, tag[1], Recvu, num, mpi_type, Talk[dir].rank_u, tag[1], Comm[DIR_ALL], &status[1]);
  }

#endif // GKC_PARALLEL_MPI
  
  return;
}


void Parallel::updateBoundaryVlasovBarrier() 
{
  // BUG what happen if we never sent a message, what does Waitall
  check(MPI_Waitall(4, Talk[DIR_X].psf_msg_req, Talk[DIR_X].msg_status) != MPI_SUCCESS ? -1 : 0, DMESG("MPI Error"));
  //if(decomposition[DIR_X] > 1) MPI_Waitall(4, Talk[DIR_X].psf_msg_req, Talk[DIR_X].msg_status);
  if(Nz > 1) MPI_Waitall(4, Talk[DIR_Z].psf_msg_req, Talk[DIR_Z].msg_status);
  if(decomposition[DIR_V] > 1) MPI_Waitall(4, Talk[DIR_V].psf_msg_req, Talk[DIR_V].msg_status);
      
  return;
}


// There should be a better way instead of defininng 2 updateNEighbours as all same the same functions
// but template arguments are different ... :(

void  Parallel::updateBoundaryFields(CComplex *SendXl, CComplex *SendXu, CComplex *RecvXl, CComplex *RecvXu, int num_X,
                                     CComplex *SendZl, CComplex *SendZu, CComplex *RecvZl, CComplex *RecvZu, int num_Z) 
{
   
  MPI_Status  msg_status [8];
  MPI_Request msg_request[8];
      
  // For X-Direction 
  MPI_Irecv(RecvXl, num_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_l, Talk[DIR_X].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[0]); 
  MPI_Isend(SendXu, num_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_u, Talk[DIR_X].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[1]);
  
  MPI_Irecv(RecvXu, num_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_u, Talk[DIR_X].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[2]); 
  MPI_Isend(SendXl, num_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_l, Talk[DIR_X].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[3]);

  // For Z-Direction (ignore if we don't use Z direction)
  if(Nz > 1) {
  MPI_Irecv(RecvZl, num_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_l, Talk[DIR_Z].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[4]); 
  MPI_Isend(SendZu, num_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_u, Talk[DIR_Z].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[5]);
      
  MPI_Irecv(RecvZu, num_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_u, Talk[DIR_Z].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[6]); 
  MPI_Isend(SendZl, num_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_l, Talk[DIR_Z].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[7]);
  }
  // Field boundries required for Vlasov solver, thus wait
  MPI_Waitall(Nz > 1 ? 8 : 4, msg_request, msg_status);
      
  return;
}


MPI_Op Parallel::getMPIOp(Op op) 
{
  
  MPI_Op mOp = MPI_OP_NULL;

  switch(op) {

    case(Op::sum ) : mOp = MPI_SUM ; break;
    case(Op::max ) : mOp = MPI_MAX ; break;
    case(Op::min ) : mOp = MPI_MIN ; break;
    case(Op::bor ) : mOp = MPI_BOR ; break;
    case(Op::band) : mOp = MPI_BAND; break;
    case(Op::lor ) : mOp = MPI_LOR ; break;
    case(Op::land) : mOp = MPI_LAND; break;

    default        : check(-1, DMESG("No such MPI operation defined"));
  }
  
  return mOp;
}

MPI_Datatype Parallel::getMPIDataType(const std::type_info &T) 
{

  MPI_Datatype type=0;

  if     (T == typeid( Complex ) ) type = MPI_DOUBLE_COMPLEX;
  else if(T == typeid(CComplex ) ) type = MPI_DOUBLE_COMPLEX;
  else if(T == typeid(double   ) ) type = MPI_DOUBLE;
  else if(T == typeid(int      ) ) type = MPI_INT;
  else if(T == typeid(long long) ) type = MPI_LONG_LONG;
  //else if(T == typeid(bool     ) ) type = MPI::BOOL;   
  else check(-1, DMESG("Such type is not defined"));
    
  return type;

}


void Parallel::getAutoDecomposition(int numCPU) 
{

  if (numCPU == 1) decomposition[:] = 1;
  else check(-1, DMESG("Not implemented"));

  return;
}
  

void Parallel::checkValidDecomposition(Setup *setup) 
{

  // Check basic decomposition sizes
  if( decomposition[DIR_X] > setup->get("Grid.Nx", 1)) check(-1, DMESG("Decomposition in x bigger than Nx"));
  if( decomposition[DIR_Z] > setup->get("Grid.Nz", 1)) check(-1, DMESG("Decomposition in z bigger than Nz"));
  if( decomposition[DIR_V] > setup->get("Grid.Nv", 1)) check(-1, DMESG("Decomposition in v bigger than Nv"));
  if( decomposition[DIR_M] > setup->get("Grid.Nm", 1)) check(-1, DMESG("Decomposition in m bigger than Nm"));
  if( decomposition[DIR_S] > setup->get("Grid.Ns", 1)) check(-1, DMESG("Decomposition in s bigger than Ns"));
  
  // check if OpenMP threads are equal decompositon
  if( decomposition[DIR_Y] != numThreads             ) check(-1, DMESG("Failed to set threads. Decomposition[Y] != numThreads. Check OMP_NUM_THREADS"));
   
  // Simple Check if reasonable values are provided for decomposition (only MPI proceeses)
  const int pNs = setup->get("Grid.Ns", 1 );
  //if((product(decomposition) != numThreads * numProcesses) && (myRank == 0)) check(-1, DMESG("Decomposition and number of processors are not equal"));
  if((__sec_reduce_mul(decomposition[DIR_X:6]) != numThreads * numProcesses) && (myRank == 0)) check(-1, DMESG("Decomposition and number of processors are not equal"));
  if(((pNs %   decomposition[DIR_S]) != 0    ) && (myRank == 0)) check(-1, DMESG("Decomposition in s have to be modulo of the total number"));

   return;
}


void Parallel::print(std::string message)
{

  if(myRank == 0) std::cout << message << std::endl;
}

void Parallel::printOn(std::ostream &output) const 
{

  output << "Parallel   | ";
      
  if      (useOpenMP) output << " OpenMP (Threads) : " << Setup::num2str(numThreads) ;
  if      (useMPI   ) {
      
  output << " MPI (processes) : " << Setup::num2str(numProcesses) << std::endl;
    
  output << "           |  Decompostion : " 
         << "X(" << decomposition[DIR_X] << ")  Y(" << decomposition[DIR_Y] << ") "
         << "Z(" << decomposition[DIR_Z] << ")  V(" << decomposition[DIR_V] << ") "
         << "M(" << decomposition[DIR_M] << ")  S(" << decomposition[DIR_S] << ") " << std::endl;
   
  } else output << std::endl;
} 


void Parallel::initData(Setup *setup, FileIO *fileIO) 
{

  hid_t parallelGroup = fileIO->newGroup("/Parallel");
   
  check(H5LTset_attribute_int(parallelGroup, ".", "X",  &decomposition[DIR_X], 1), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_int(parallelGroup, ".", "Y",  &decomposition[DIR_Y], 1), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_int(parallelGroup, ".", "Z",  &decomposition[DIR_Z], 1), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_int(parallelGroup, ".", "V",  &decomposition[DIR_V], 1), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_int(parallelGroup, ".", "M",  &decomposition[DIR_M], 1), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_int(parallelGroup, ".", "S",  &decomposition[DIR_S], 1), DMESG("H5LTset_attribute"));
  
  // save MPI library information
  int version, subversion;
  MPI_Get_version(&version, &subversion);
  
  H5LTset_attribute_string(parallelGroup, ".", "MPI Version",  (Setup::num2str(version) + "." + Setup::num2str(subversion)).c_str());
  
  H5Gclose(parallelGroup);
}


int Parallel::getNumberOfWorkers(int dir) 
{
  int numWorkers = 0;
  MPI_Comm_size(Comm[dir],&numWorkers);
  return numWorkers;
}
    
   
int Parallel::getWorkerID(int dir) 
{
  int rankWorker = 0;
  MPI_Comm_rank(Comm[dir],&rankWorker);
  return rankWorker;
}

void Parallel::barrier(int dir) 
{
  MPI_Barrier(Comm[dir]);
}


void Parallel::printProcessID()
{
  std::cout << " Parallel rank  " << myRank 
            << " Thread   id    " << threadID
  << " Decomposition  "  
         << "X(" << decomposition[DIR_X] << ")  Y(" << decomposition[DIR_Y] << ") "
         << "Z(" << decomposition[DIR_Z] << ")  V(" << decomposition[DIR_V] << ") "
         << "M(" << decomposition[DIR_M] << ")  S(" << decomposition[DIR_S] << ") " 
         << std::endl << std::flush;


}

void Parallel::setThreadID()
{
#ifndef GKC_PARALLEL_OPENMP
  threadID = 0;
#else 
  threadID = omp_get_thread_num();
#endif    
}
