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

// global process rank needed for output
extern int process_rank;

// MPI_ERror handler so we can use check and backtrace
void check_mpi(MPI_Comm *comm, int *err_code, ...) {

  char string[MPI_MAX_ERROR_STRING]; int len;
  MPI_Error_string(*err_code, string, &len);

  if (*err_code != MPI_SUCCESS) check(-1, DMESG(std::string(string)), true);

}



Parallel::Parallel(Setup *setup) 
{
  // initialize some basic parameter
  myRank = 0;  numThreads = 1; numProcesses = 1, master_rank = 0; 
  useOpenMP = false; useMPI = true;
  Coord[:] = 0;

  decomposition[:] = 1;

#ifdef GKC_PARALLEL_OPENMP
  useOpenMP = true;
  // if OpenMP is enabled, we decompose over y direction
  #pragma omp parallel
  {
    numThreads = omp_get_num_threads();
  }
#endif // GKC_PARALLEL_OPENMP

   for(int d=DIR_X;d<DIR_SIZE;d++) Comm[d] = MPI_COMM_NULL;

   //////////////////////// Set Message tags, enumerate through ////////////////////////
   int i = 14665; // some random number to not to interfere with other
   // do its in one loop ?
   for(int dir = DIR_X; dir <= DIR_S; dir++) {
        Talk[dir].psf_msg_tag[0] = ++i;
        Talk[dir].psf_msg_tag[1] = ++i;
   }
   for(int dir = DIR_X; dir <= DIR_S; dir++) {
        Talk[dir].phi_msg_tag[0] = ++i;
        Talk[dir].phi_msg_tag[1] = ++i;
   }
   
   // MPI-2 standard allows to pass NULL for (&argc, &argv)
   // valgrind complains about unilizaedMPI_Init     (&setup->argc, &setup->argv);
   MPI_Init     (NULL, NULL);

   MPI_Comm_size(MPI_COMM_WORLD, &numProcesses); 
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // set automatic decomposition
  std::vector<std::string> decomp = Setup::split(setup->get("Parallel.Decomposition","Auto"), ":");
  // check also if decomp is not longer than N<=6
  if(decomp[0] == "Auto")  getAutoDecomposition(numProcesses);
  else for(unsigned int dir=DIR_X; dir < decomp.size() && (dir <= DIR_S); dir++)  decomposition[dir] =  atoi(decomp[dir].c_str());
 
  // Note : if OpenMP is enabled OpenMP, we decompose in Y in OpenMP threads (not clean solution tough)
#ifdef GKC_PARALLEL_OPENMP
  #pragma omp parallel
  {
      omp_set_num_threads(numThreads);
  } 
#endif

   checkValidDecomposition(setup);
   
   // need to set to 1 for MPI
   int dec_Y = decomposition[DIR_Y]; decomposition[DIR_Y] = 1; 
   
    
   /////////////////// Create 6-D Cartesian Grid ////////////////////
   int periods[6] = { 1, 1, 1, 0, 0, 0};
   MPI_Cart_create(MPI_COMM_WORLD, 6, decomposition, periods, 1, &Comm[DIR_ALL]);

   // Get Cart coordinates to set 
   MPI_Comm_rank(Comm[DIR_ALL],&myRank); 
   process_rank = myRank; // process rank is global, needed ?

   MPI_Cart_coords(Comm[DIR_ALL], myRank, 6, Coord);
 
   // Get combined positions 
   Coord[DIR_VMS]   = ((Coord[DIR_V] == 0) && (Coord[DIR_M] == 0) && (Coord[DIR_S] == 0)) ? 0 : -1;
   Coord[DIR_VM ]   = ((Coord[DIR_V] == 0) && (Coord[DIR_M] == 0))                        ? 0 : -1;
   Coord[DIR_MS ]   = ((Coord[DIR_M] == 0) && (Coord[DIR_S] == 0))                        ? 0 : -1;

    ///////////////////  Set SubCommunicators for various directions //////////////////////////
    // Note , we do not need to translate ranks, because dirMaster is always used with appropriate Comm
    
    int coord_master[6] = { 0, 0, 0, 0, 0, 0 };
 
    // Communicator for X
    int remain_dims_X[6] = { true, false, false, false, false, false };         
    MPI_Cart_sub(Comm[DIR_ALL], remain_dims_X, &Comm[DIR_X]);
    MPI_Cart_rank  (Comm[DIR_X], coord_master, &dirMaster[DIR_X]);
    
    // Communicator for Y
    int remain_dims_Y[6] = { false, true, false, false, false, false };         
    MPI_Cart_sub(Comm[DIR_ALL], remain_dims_Y, &Comm[DIR_Y]);
    MPI_Cart_rank  (Comm[DIR_Y], coord_master, &dirMaster[DIR_Y]);
    
    // Communicator for Z
    int remain_dims_Z[6] = { false, false, true, false, false, false };         
    MPI_Cart_sub(Comm[DIR_ALL], remain_dims_Z, &Comm[DIR_Z]);
    MPI_Cart_rank  (Comm[DIR_Z], coord_master, &dirMaster[DIR_Z]);
    
    // Communicator for V
    int remain_dim_V[6]   = {false, false, false, true, false ,false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_V, &Comm[DIR_V]);
    MPI_Cart_rank  (Comm[DIR_V], coord_master, &dirMaster[DIR_V]);

    // Communicator for M
    int remain_dim_M[6] = {false, false, false, false, true ,false};
    MPI_Cart_sub  (Comm[DIR_ALL], remain_dim_M, &Comm[DIR_M]); 
    MPI_Cart_rank (Comm[DIR_M]  , coord_master, &dirMaster[DIR_M]);

    // Communicator for S
    int remain_dim_S[6] = {false, false, false, false, false, true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_S, &Comm[DIR_S]); 
    MPI_Cart_rank  (Comm[DIR_S], coord_master, &dirMaster[DIR_S]);
   
    // Communicator for M&S (used by gyro-averaging solver)
    int remain_dim_MS[6] = {false, false, false, false, true, true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_MS, &Comm[DIR_MS]); 
    MPI_Cart_rank  (Comm[DIR_MS], coord_master, &dirMaster[DIR_MS]);
   
    // For real space communication
    int remain_dim_XYZ[6] = {true, true, true, false, false, false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_XYZ, &Comm[DIR_XYZ]); 
    MPI_Cart_rank  (Comm[DIR_XYZ], coord_master, &dirMaster[DIR_XYZ]);
    
    int remain_dim_VMS[6] = {false, false, false, true, true ,true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_VMS, &Comm[DIR_VMS]);
    MPI_Cart_rank  (Comm[DIR_VMS], coord_master, &dirMaster[DIR_VMS]);
   
    // Velocity space communiocation
    int remain_dim_VM[6] = { false, false, false, true, true, false };
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_VM, &Comm[DIR_VM]);
    MPI_Cart_rank  (Comm[DIR_VM], coord_master, &dirMaster[DIR_VM]);
    
    // Phase-space of individual species
    int remain_dim_XYZVM[6] = { true, true, true, true, true ,false };
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_XYZVM, &Comm[DIR_XYZVM]);
    MPI_Cart_rank  (Comm[DIR_XYZVM], coord_master, &dirMaster[DIR_XYZVM]);
    
    ////// Get ranks of neigbouring processes for Vlasov equation boundary exchange //////

    int rank = myRank; // is this necessary ?!
    
    // setup communcation ranks for X
    MPI_Cart_shift(Comm[DIR_ALL], DIR_X, 1, &rank, &Talk[DIR_X].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], DIR_X,-1, &rank, &Talk[DIR_X].rank_l);
  
    // setup communcation ranks for Y
    MPI_Cart_shift(Comm[DIR_ALL], DIR_Y, 1, &rank, &Talk[DIR_Y].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], DIR_Y,-1, &rank, &Talk[DIR_Y].rank_l);
 
    // setup communcation ranks for Z
    MPI_Cart_shift(Comm[DIR_ALL], DIR_Z, 1, &rank, &Talk[DIR_Z].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], DIR_Z,-1, &rank, &Talk[DIR_Z].rank_l);

    // V (not periodic), returns NULL_COMM at the domain ends
    MPI_Cart_shift(Comm[DIR_ALL], DIR_V, 1, &rank, &Talk[DIR_V].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], DIR_V,-1, &rank, &Talk[DIR_V].rank_l);
  
    // S (not periodic), returns NULL_COMM at the domain ends
    MPI_Cart_shift(Comm[DIR_ALL], DIR_M, 1, &rank, &Talk[DIR_M].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], DIR_M,-1, &rank, &Talk[DIR_M].rank_l);
  
    // M (not periodic), returns NULL_COMM at the domain ends
    MPI_Cart_shift(Comm[DIR_ALL], DIR_S, 1, &rank, &Talk[DIR_S].rank_u);
    MPI_Cart_shift(Comm[DIR_ALL], DIR_S,-1, &rank, &Talk[DIR_S].rank_l);
 
    /////////////////// Setup own error handle,  this needed to enable backtracing when debugging
    MPI_Errhandler my_errhandler;
    MPI_Errhandler_create(&check_mpi, &my_errhandler);
    MPI_Errhandler_set(Comm[DIR_ALL], my_errhandler);

}



Parallel::~Parallel() 
{
   
     // free MPI communicators (Comm[] was initialized with MPI_COMM_NULL)
      for(int n=DIR_X; n<DIR_SIZE; n++) if(Comm[n] != MPI_COMM_NULL) MPI_Comm_free(&Comm[n]);

     check(MPI_Finalize(), DMESG("MPI_Finalize"));
 
     return;
}

void Parallel::updateNeighboursBarrier() 
{

    // BUG what happen if we never sent a message, what does Waitall
     
    //if(decomposition[DIR_X] > 1) MPI_Waitall(4, Talk[DIR_X].psf_msg_request, Talk[DIR_X].msg_status);
      MPI_Waitall(4, Talk[DIR_X].psf_msg_request, Talk[DIR_X].msg_status);
    if(decomposition[DIR_Y] > 1) MPI_Waitall(4, Talk[DIR_Y].psf_msg_request, Talk[DIR_Y].msg_status);
    if(decomposition[DIR_Z] > 1) MPI_Waitall(4, Talk[DIR_Z].psf_msg_request, Talk[DIR_Z].msg_status);
    if(decomposition[DIR_V] > 1) MPI_Waitall(4, Talk[DIR_V].psf_msg_request, Talk[DIR_V].msg_status);
  //if(decomposition[DIR_M] > 1) MPI_Waitall(4, Talk[DIR_M].psf_msg_request, Talk[DIR_M].msg_status);
  //if(decomposition[DIR_S] > 1) MPI_Waitall(4, Talk[DIR_S].psf_msg_request, Talk[DIR_S].msg_status);
      
    return;
}


// There should be a better way instead of defininng 2 updateNEighbours as all same the same functions
// but template arguments are different ... :(

void  Parallel::updateBoundaryFields(CComplex *SendXl, CComplex *SendXu, CComplex *RecvXl, CComplex *RecvXu, int numElements_X,
                                     CComplex *SendZl, CComplex *SendZu, CComplex *RecvZl, CComplex *RecvZu, int numElements_Z) 
//void Parallel::updateNeighbours(Array6C  SendXl, Array6C  SendXu, Array6C  SendYl, Array6C  SendYu, Array6C SendZl, Array6C SendZu, 
//                                Array6C  RecvXl, Array6C  RecvXu, Array6C  RecvYl, Array6C  RecvYu, Array6C RecvZl, Array6C RecvZu) 
{
      MPI_Status  msg_status[8];
      MPI_Request msg_request[8];
      
      // For X-Direction 
      MPI_Irecv(RecvXl, numElements_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_l, Talk[DIR_X].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[0]); 
      MPI_Isend(SendXu, numElements_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_u, Talk[DIR_X].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[1]);
      
      MPI_Irecv(RecvXu, numElements_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_u, Talk[DIR_X].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[2]); 
      MPI_Isend(SendXl, numElements_X, MPI_DOUBLE_COMPLEX, Talk[DIR_X].rank_l, Talk[DIR_X].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[3]);

      // For Y-Direction
//      MPI_Irecv(RecvYl, RecvYl.numElements(), MPI_DOUBLE_COMPLEX, Talk[DIR_Y].rank_l, Talk[DIR_Y].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[4]); 
//      MPI_Isend(SendYu, SendYu.numElements(), MPI_DOUBLE_COMPLEX, Talk[DIR_Y].rank_u, Talk[DIR_Y].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[5]);
      
//      MPI_Irecv(RecvYu, RecvYu.numElements(), MPI_DOUBLE_COMPLEX, Talk[DIR_Y].rank_u, Talk[DIR_Y].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[6]); 
//      MPI_Isend(SendYl, SendYl.numElements(), MPI_DOUBLE_COMPLEX, Talk[DIR_Y].rank_l, Talk[DIR_Y].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[7]);
      
      // For Z-Direction
      MPI_Irecv(RecvZl, numElements_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_l, Talk[DIR_Z].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[4]); 
      MPI_Isend(SendZu, numElements_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_u, Talk[DIR_Z].phi_msg_tag[0], Comm[DIR_ALL], &msg_request[5]);
      
      MPI_Irecv(RecvZu, numElements_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_u, Talk[DIR_Z].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[6]); 
      MPI_Isend(SendZl, numElements_Z, MPI_DOUBLE_COMPLEX, Talk[DIR_Z].rank_l, Talk[DIR_Z].phi_msg_tag[1], Comm[DIR_ALL], &msg_request[7]);

      // Ok let's wait here ....
      MPI_Waitall(8, msg_request, msg_status);
      
      return;
}

MPI_Op Parallel::getMPIOp(Op op) 
{

    MPI_Op mOp = MPI_OP_NULL;

    switch(op) {
      case(Op::SUM)  : mOp = MPI_SUM ; break;
      case(Op::MAX)  : mOp = MPI_MAX ; break;
      case(Op::MIN)  : mOp = MPI_MIN ; break;
      case(Op::BOR ) : mOp = MPI_BOR ; break;
      case(Op::BAND) : mOp = MPI_BAND; break;
      case(Op::LOR ) : mOp = MPI_LOR ; break;
      case(Op::LAND) : mOp = MPI_LAND; break;

        default       : check(-1, DMESG("No such MPI operation defined"));
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
    else if(T == typeid(bool     ) ) type = MPI::BOOL;   
    else check(-1, DMESG("Such type is not defined"));
    
   return type;
}


void Parallel::getAutoDecomposition(int numCPU) 
{

    if (numCPU == 1) decomposition[:] = 1;
    else check(-1, DMESG("Not implemented"));

    return;

};
  

bool Parallel::checkValidDecomposition(Setup *setup) 
{

   // Check basic decomposition sizes
   if( decomposition[DIR_X] > setup->get("Grid.Nx", 1)) check(-1, DMESG("Decomposition in x bigger than Nx"));
// no need to check y-decomposition (OpenMP parallelization)   
// if( decomposition(DIR_Y) > setup->get("Grid.Ny", 1)) check(-1, DMESG("Decomposition in y bigger than Ny"));
   if( decomposition[DIR_Z] > setup->get("Grid.Nz", 1)) check(-1, DMESG("Decomposition in z bigger than Nz"));
   if( decomposition[DIR_V] > setup->get("Grid.Nv", 1)) check(-1, DMESG("Decomposition in v bigger than Nv"));
   if( decomposition[DIR_M] > setup->get("Grid.Nm", 1)) check(-1, DMESG("Decomposition in m bigger than Nm"));
   if( decomposition[DIR_S] > setup->get("Grid.Ns", 1)) check(-1, DMESG("Decomposition in s bigger than Ns"));
   
   // Simple Check if reasonable values are provided for decomposition (only MPI proceeses)
   const int pNs = setup->get("Grid.Ns", 1 );
   //if((product(decomposition) != numThreads * numProcesses) && (myRank == 0)) check(-1, DMESG("Decomposition and number of processors are not equal"));
   if((__sec_reduce_mul(decomposition[DIR_X:6]) != numThreads * numProcesses) && (myRank == 0)) check(-1, DMESG("Decomposition and number of processors are not equal"));
   if(((pNs %   decomposition[DIR_S]) != 0    ) && (myRank == 0)) check(-1, DMESG("Decomposition in s have to be modulo of the total number"));


    return GKC_SUCCESS;
};


void Parallel::print(std::string message)
{

    if(myRank == 0) std::cout << message << std::endl;

}


void Parallel::printOn(std::ostream &output) const {

       output << "Parallel   | ";
      
      if      (useOpenMP) output << " OpenMP (Threads) : " << Setup::num2str(numThreads) ;
      if      (useMPI   ) {
      
           if (decomposition[DIR_X] == 0) output <<  "Automatic" << std::endl;
           else  output <<  " MPI (processes) : " << Setup::num2str(numProcesses) << std::endl;

       
      output << "           |  Decompostion : " 
        << "X(" << decomposition[DIR_X] << ")  Y(" << decomposition[DIR_Y] << ") "
        << "Z(" << decomposition[DIR_Z] << ")  V(" << decomposition[DIR_V] << ") "
        << "M(" << decomposition[DIR_M] << ")  S(" << decomposition[DIR_S] << ") " << std::endl;
   
       } else output << std::endl;
} 

int Parallel::getNumberOfWorkers(int dir) 
{
        int numWorkers = 0;
        MPI_Comm_size(Comm[dir],&numWorkers);
        return numWorkers;
   
};
    
   
int Parallel::getWorkerID(int dir) 
{
        int rankWorker = 0;
        MPI_Comm_rank(Comm[dir],&rankWorker);
        return rankWorker;

};

void Parallel::barrier(int dir) 
{
        MPI_Barrier(Comm[dir]);
};
