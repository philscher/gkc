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


#ifdef GKC_PARALLEL_MPI
// MPI_ERror handler so we can use check and backtrace
void check_mpi(MPI_Comm *comm, int *err_code, ...) {

  char string[MPI_MAX_ERROR_STRING]; int len;
  MPI_Error_string(*err_code, string, &len);

  if (*err_code != MPI_SUCCESS) check(-1, DMESG(std::string(string)), true);

}
#endif // GKC_PARALLEL_MPI



Parallel::Parallel(Setup *setup) 
{
  // initialize some basic parameter
  myRank = 0;  numThreads = 1; numProcesses = 1, master_rank = 0; numGPUs = 0; 
  useOpenMP = false; useMPI = false; useOpenCL = false; 
  Coord.resize(Range(DIR_X,DIR_SIZE-1)); Coord(Range(DIR_X, DIR_SIZE-1)) = 0; //Coord(Range(DIR_S+1, DIR_SIZE-1)) = -1;
  master_process_id = setup->get("Helios.Process_ID", 0);

#ifdef PARALLEL_OPENMP

      useOpenMP = true;

      // if OpenMP is enabled, we decompose over y direction

      #pragma omp parallel
      {
        numThreads = omp_get_num_threads();
      }
#endif

    /////////////////////////////////////////////////////////////////////////////////////////
    // Check domain decomposition
    decomposition.resize(Range(DIR_X,DIR_S)); decomposition = 1;
     /////////////////////////////////////////////////////////////////////////////////////////////

#ifdef GKC_PARALLEL_MPI
   useMPI = true;
   for(int d=DIR_X;d<DIR_SIZE;d++) Comm[d] = MPI_COMM_NULL;
   Talk.resize(Range(DIR_X,DIR_S));

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // Set Message tags, enumerate through
   int i = 0;
   
   for(int dir = DIR_X; dir <= DIR_S; dir++) {
        Talk(dir).psf_msg_tag[0] = ++i;
        Talk(dir).psf_msg_tag[1] = ++i;
   }
   for(int dir = DIR_X; dir <= DIR_S; dir++) {
        Talk(dir).phi_msg_tag[0] = ++i;
        Talk(dir).phi_msg_tag[1] = ++i;
   }
   
   MPI_Init(&setup->argc,&setup->argv); 
   MPI_Comm_size(MPI_COMM_WORLD,&numProcesses); 
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

  // set automatic decomposition
  std::vector<std::string> decomp = Setup::split(setup->get("Parallel.Decomposition","Auto"), ":");
  if(decomp[0] == "Auto")                              decomposition      = getAutoDecomposition(numProcesses);
  else for(unsigned int dir=DIR_X; dir < decomp.size(); dir++)  decomposition(dir) =  atoi(decomp[dir].c_str());
 
  // Not if OpenMP is enabled OpenMP, we decompose in Y in OpenMP threads (not clean solution tough)
#ifdef PARALLEL_OPENMP
      numThreads = decomposition(DIR_Y); 
      #pragma omp parallel
      {
        omp_set_num_threads(numThreads);
      } 
      // need to set to 1 for MPI
      decomposition(DIR_Y) = 1; 
#endif



  checkValidDecomposition(setup, decomposition);

   
   // Create 6-D Cartesian Grid
   int periods[6] = { 1, 1, 1, 0, 0, 0};
   MPI_Cart_create(MPI_COMM_WORLD, 6, decomposition.data(), periods, 1, &Comm[DIR_ALL]);

   // we let even fatal errors return, but we check the return code of each send
   // this needed to enable backtracing when debugging
   MPI_Errhandler my_errhandler;
   MPI_Errhandler_create(&check_mpi, &my_errhandler);
   MPI_Errhandler_set(Comm[DIR_ALL], my_errhandler);
   // Get Cart coordinates to set 
   MPI_Comm_rank(Comm[DIR_ALL],&myRank); process_rank = myRank;
   MPI_Cart_coords(Comm[DIR_ALL], myRank, 6, Coord.data());
 

   Coord(DIR_VMS)   = ((Coord(DIR_V) == 0) && (Coord(DIR_M) == 0) && (Coord(DIR_S) == 0)) ? 0 : -1;
   Coord(DIR_VM )   = ((Coord(DIR_V) == 0) && (Coord(DIR_M) == 0))                        ? 0 : -1;
   Coord(DIR_MS )   = ((Coord(DIR_M) == 0) && (Coord(DIR_S) == 0))                        ? 0 : -1;
   Coord(DIR_YZVMS) = sum(Coord(Range(DIR_Y, DIR_S))        == 0)                         ? 0 : -1;
   Coord(DIR_XYZVM) = sum(Coord(Range(DIR_X, DIR_M))        == 0)                         ? 0 : -1;

   ////////////////////////////////////////////////////////////////////////////////////////

//   if     (FFT_DECOMP == DECOMP_NO ) Coord(DIR_FFT) = ((sum(Coord(Range(DIR_X,DIR_Y)))+Coord(DIR_V)) == 0) ? 0 : -1;
 //  else if(FFT_DECOMP == DECOMP_X  ) Coord(DIR_FFT) = ((sum(Coord(Range(DIR_Y,DIR_Y)))+Coord(DIR_V)) == 0) ? 0 : -1;
 //  else   check(-1, DMESG("Spatial Decomposition strategy for FFT does not exist"));

    //  ******************  Set SubCommunicators for (V,M,S & Phi) *********************** //
    
    int coord_master[6] = { 0, 0, 0, 0, 0, 0 };
 
    // we do not need to translate ranks, because dirMaster is always used  with appropriate Comm
  
    /* 
   MPI_Group CommV_group, CommALL_group;
   MPI_Comm_group(Comm[DIR_ALL],   &CommALL_group);
   MPI_Comm_group(Comm[DIR_V  ],  &CommV_group);
*/   

    int remain_dims_X[6] = { true, false, false, false, false, false };         
    MPI_Cart_sub(Comm[DIR_ALL], remain_dims_X, &Comm[DIR_X]);
    MPI_Cart_rank  (Comm[DIR_X], coord_master, &dirMaster[DIR_X]);
    
    int remain_dims_Y[6] = { false, true, false, false, false, false };         
    MPI_Cart_sub(Comm[DIR_ALL], remain_dims_Y, &Comm[DIR_Y]);
    MPI_Cart_rank  (Comm[DIR_Y], coord_master, &dirMaster[DIR_Y]);
    
    int remain_dim_V[6]   = {false, false, false, true, false ,false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_V, &Comm[DIR_V]);
    MPI_Cart_rank  (Comm[DIR_V], coord_master, &dirMaster[DIR_V]);


    int remain_dim_M[6] = {false, false, false, false, true ,false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_M, &Comm[DIR_M]); 
    MPI_Cart_rank  (Comm[DIR_M], coord_master, &dirMaster[DIR_M]);


    int remain_dim_S[6] = {false, false, false, false, false, true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_S, &Comm[DIR_S]); 
    MPI_Cart_rank  (Comm[DIR_S], coord_master, &dirMaster[DIR_S]);
    
    int remain_dim_MS[6] = {false, false, false, false, true, true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_MS, &Comm[DIR_MS]); 
    MPI_Cart_rank  (Comm[DIR_MS], coord_master, &dirMaster[DIR_MS]);
    
    int remain_dim_XYZ[6] = {true, true, true, false, false, false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_XYZ, &Comm[DIR_XYZ]); 
    MPI_Cart_rank  (Comm[DIR_XYZ], coord_master, &dirMaster[DIR_XYZ]);
    
    int remain_dim_YZ[6] = {false, true, true, false, false, false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_YZ, &Comm[DIR_YZ]); 
    MPI_Cart_rank  (Comm[DIR_YZ], coord_master, &dirMaster[DIR_YZ]);
    
    int remain_dim_VMS[6] = {false, false, false, true, true ,true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_VMS, &Comm[DIR_VMS]);
    MPI_Cart_rank  (Comm[DIR_VMS], coord_master, &dirMaster[DIR_VMS]);
    
    int remain_dim_VM[6] = {false, false, false, true, true, false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_VM, &Comm[DIR_VM]);
    MPI_Cart_rank  (Comm[DIR_VM], coord_master, &dirMaster[DIR_VM]);
    
    int remain_dim_YZVMS[6] = {false, true, true, true, true ,true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_YZVMS, &Comm[DIR_YZVMS]);
    MPI_Cart_rank  (Comm[DIR_YZVMS], coord_master, &dirMaster[DIR_YZVMS]);
    
    int remain_dim_XYZVM[6] = {true, true, true, true, true ,false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_XYZVM, &Comm[DIR_XYZVM]);
    MPI_Cart_rank  (Comm[DIR_XYZVM], coord_master, &dirMaster[DIR_XYZVM]);
    
    int remain_dim_XMS[6] = {true, false, false, false, true , true};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_XMS, &Comm[DIR_XMS]);
    MPI_Cart_rank  (Comm[DIR_XMS], coord_master, &dirMaster[DIR_XMS]);
    
    int remain_dim_XM[6] = {true, false, false, false, true ,false};
    MPI_Cart_sub(Comm[DIR_ALL], remain_dim_XM, &Comm[DIR_XM]);
    MPI_Cart_rank  (Comm[DIR_XM], coord_master, &dirMaster[DIR_XM]);


   // because ranks changed withing groups, we need to translate them between the communicators
//   MPI_Group CommCart_group, CommSlaves_group;
//   MPI_Comm_group(CommSlaves,   &CommSlaves_group);
//   MPI_Comm_group(Comm[DIR_ALL],  &CommCart_group);
   
//   MPI_Group_translate_ranks(CommCart_group, 1, &rankFFTMaster,  CommSlaves_group, &rankFFTMaster);


  // Get Communicators to our positions, use MPI_CART_SHIFT for that !
  int rank = myRank;
  MPI_Cart_shift(Comm[DIR_ALL], DIR_X, 1, &rank, &Talk(DIR_X).rank_u);
  MPI_Cart_shift(Comm[DIR_ALL], DIR_X,-1, &rank, &Talk(DIR_X).rank_l);
  
  // setup communcation ranks for Y
  MPI_Cart_shift(Comm[DIR_ALL], DIR_Y, 1, &rank, &Talk(DIR_Y).rank_u);
  MPI_Cart_shift(Comm[DIR_ALL], DIR_Y,-1, &rank, &Talk(DIR_Y).rank_l);
 
  // setup communcation ranks for Z
  MPI_Cart_shift(Comm[DIR_ALL], DIR_Z, 1, &rank, &Talk(DIR_Z).rank_u);
  MPI_Cart_shift(Comm[DIR_ALL], DIR_Z,-1, &rank, &Talk(DIR_Z).rank_l);

  // V (not periodic), returns NULL_COMM at the domain ends
  MPI_Cart_shift(Comm[DIR_ALL], DIR_V, 1, &rank, &Talk(DIR_V).rank_u);
  MPI_Cart_shift(Comm[DIR_ALL], DIR_V,-1, &rank, &Talk(DIR_V).rank_l);
  
  // S (not periodic), returns NULL_COMM at the domain ends
  MPI_Cart_shift(Comm[DIR_ALL], DIR_M, 1, &rank, &Talk(DIR_M).rank_u);
  MPI_Cart_shift(Comm[DIR_ALL], DIR_M,-1, &rank, &Talk(DIR_M).rank_l);
  
  // M (not periodic), returns NULL_COMM at the domain ends
  MPI_Cart_shift(Comm[DIR_ALL], DIR_S, 1, &rank, &Talk(DIR_S).rank_u);
  MPI_Cart_shift(Comm[DIR_ALL], DIR_S,-1, &rank, &Talk(DIR_S).rank_l);
 
  // distributed unified id to all processes
  master_process_id = (int) collect((double) (myRank == 0) ? master_process_id : 0, OP_SUM, DIR_ALL);

#endif // GKC_PARALLEL_MPI


}



Parallel::~Parallel() {

#ifdef GKC_PARALLEL_MPI 


     

     // free supcommunicators
     //for(int d=DIR_X;d<DIR_SIZE;d++) Comm[d] = MPI_COMM_NULL;
     check(MPI_Comm_free(&Comm[DIR_VMS]   ), DMESG("MPI_Comm_free(Comm_S)"   ));
     check(MPI_Comm_free(&Comm[DIR_MS ]   ), DMESG("MPI_Comm_free(Comm_S)"   ));
     check(MPI_Comm_free(&Comm[DIR_V  ]   ), DMESG("MPI_Comm_free(Comm_S)"   ));
     check(MPI_Comm_free(&Comm[DIR_XYZ] ), DMESG("MPI_Comm_free(Comm_XYZ)" ));
     check(MPI_Comm_free(&Comm[DIR_ALL] ), DMESG("MPI_Comm_free(Comm_XYZ)" ));

     check(MPI_Finalize(), DMESG("MPI_Finalize"));
#endif // GKC_PARALLEL_MPI
 }

/* 

// There should be a better way instead of defininng 2 updateNEighbours as all same the same functions
// but teplate arguments are different ... :(, if the border is reached we send recv value to zero !
int Parallel::updateNeighbours(Array6z  Sendu, Array6z  Sendl, Array6z  Recvu, Array6z  Recvl, int dir, bool nonBlocking)
{
#ifdef GKC_PARALLEL_MPI
    if(nonBlocking) { 
        MPI_Irecv(Recvl.data(), Recvl.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_l, Talk(dir).psf_msg_tag[0], Comm[DIR_ALL], &Talk(dir).psf_msg_request[1]);
        MPI_Isend(Sendu.data(), Sendu.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_u, Talk(dir).psf_msg_tag[0], Comm[DIR_ALL], &Talk(dir).psf_msg_request[0]);
        if(Talk(dir).rank_u == MPI_PROC_NULL)   Recvl = 0.e0;
     
        MPI_Irecv(Recvu.data(), Recvu.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_u, Talk(dir).psf_msg_tag[1], Comm[DIR_ALL], &Talk(dir).psf_msg_request[3]); 
        MPI_Isend(Sendl.data(), Sendl.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_l, Talk(dir).psf_msg_tag[1], Comm[DIR_ALL], &Talk(dir).psf_msg_request[2]);
        if(Talk(dir).rank_l == MPI_PROC_NULL)   Recvu = 0.e0;
    }
    else {
        MPI_Sendrecv(Sendu.data(), Sendu.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_u, Talk(dir).psf_msg_tag[1], 
                     Recvl.data(), Recvl.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_l, Talk(dir).psf_msg_tag[0], Comm[dir], Talk(dir).msg_status);
        MPI_Sendrecv(Sendl.data(), Sendl.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_l, Talk(dir).psf_msg_tag[1], 
                     Recvu.data(), Recvu.numElements(), MPI_DOUBLE_COMPLEX, Talk(dir).rank_u, Talk(dir).psf_msg_tag[0], Comm[dir], Talk(dir).msg_status);
    }
      return GKC_SUCCESS;
}
 */

int Parallel::updateNeighboursBarrier() {
    // BUG what happen if we never sent a message, what does Waitall
#ifdef GKC_PARALLEL_MPI 
     if(decomposition(DIR_X) > 1) MPI_Waitall(4, Talk(DIR_X).psf_msg_request, Talk(DIR_X).msg_status);
     if(decomposition(DIR_Y) > 1) MPI_Waitall(4, Talk(DIR_Y).psf_msg_request, Talk(DIR_Y).msg_status);
     if(decomposition(DIR_Z) > 1) MPI_Waitall(4, Talk(DIR_Z).psf_msg_request, Talk(DIR_Z).msg_status);
     if(decomposition(DIR_V) > 1) MPI_Waitall(4, Talk(DIR_V).psf_msg_request, Talk(DIR_V).msg_status);
//     if(decomposition & DECOMP_M) MPI_Waitall(4, Talk(DIR_M).psf_msg_request, Talk(DIR_M).msg_status);
//     if(decomposition & DECOMP_S) MPI_Waitall(4, Talk(DIR_S).psf_msg_request, Talk(DIR_S).msg_status);
#endif // GKC_PARALLEL_MPI
      return GKC_SUCCESS;

}


// There should be a better way instead of defininng 2 updateNEighbours as all same the same functions
// but teplate arguments are different ... :(
int Parallel::updateNeighbours(Array6z  SendXl, Array6z  SendXu, Array6z  SendYl, Array6z  SendYu, Array6z SendZl, Array6z SendZu, 
                               Array6z  RecvXl, Array6z  RecvXu, Array6z  RecvYl, Array6z  RecvYu, Array6z RecvZl, Array6z RecvZu) 
{
#ifdef GKC_PARALLEL_MPI 
      MPI_Status  msg_status[12];
      MPI_Request msg_request[12];
      
      // For X-Direction 
      MPI_Irecv(RecvXl.data(), RecvXl.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_X).rank_l, Talk(DIR_X).phi_msg_tag[0], Comm[DIR_ALL], &msg_request[5]); 
      MPI_Isend(SendXu.data(), SendXu.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_X).rank_u, Talk(DIR_X).phi_msg_tag[0], Comm[DIR_ALL], &msg_request[4]);
      
      MPI_Irecv(RecvXu.data(), RecvXu.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_X).rank_u, Talk(DIR_X).phi_msg_tag[1], Comm[DIR_ALL], &msg_request[7]); 
      MPI_Isend(SendXl.data(), SendXl.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_X).rank_l, Talk(DIR_X).phi_msg_tag[1], Comm[DIR_ALL], &msg_request[6]);

      // For Y-Direction
      MPI_Irecv(RecvYl.data(), RecvYl.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Y).rank_l, Talk(DIR_Y).phi_msg_tag[0], Comm[DIR_ALL], &msg_request[1]); 
      MPI_Isend(SendYu.data(), SendYu.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Y).rank_u, Talk(DIR_Y).phi_msg_tag[0], Comm[DIR_ALL], &msg_request[0]);
      
      MPI_Irecv(RecvYu.data(), RecvYu.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Y).rank_u, Talk(DIR_Y).phi_msg_tag[1], Comm[DIR_ALL], &msg_request[3]); 
      MPI_Isend(SendYl.data(), SendYl.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Y).rank_l, Talk(DIR_Y).phi_msg_tag[1], Comm[DIR_ALL], &msg_request[2]);
      
      // For Z-Direction
      MPI_Irecv(RecvZl.data(), RecvZl.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Z).rank_l, Talk(DIR_Z).phi_msg_tag[0], Comm[DIR_ALL], &msg_request[9]); 
      MPI_Isend(SendZu.data(), SendZu.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Z).rank_u, Talk(DIR_Z).phi_msg_tag[0], Comm[DIR_ALL], &msg_request[8]);
      
      MPI_Irecv(RecvZu.data(), RecvZu.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Z).rank_u, Talk(DIR_Z).phi_msg_tag[1], Comm[DIR_ALL], &msg_request[11]); 
      MPI_Isend(SendZl.data(), SendZl.numElements(), MPI_DOUBLE_COMPLEX, Talk(DIR_Z).rank_l, Talk(DIR_Z).phi_msg_tag[1], Comm[DIR_ALL], &msg_request[10]);


      // Ok let's wait here ....
      MPI_Waitall(12, msg_request, msg_status);
#endif // GKC_PARALLEL_MPI
      return GKC_SUCCESS;
}

#ifdef GKC_PARALLEL_MPI

MPI_Op Parallel::getMPIOp(int op) {
    MPI_Op mOp = MPI_OP_NULL;
    switch(op) {
        case(OP_SUM) : mOp = MPI_SUM; break;
        case(OP_MAX) : mOp = MPI_MAX; break;
        case(OP_MIN) : mOp = MPI_MIN; break;
        case(OP_BOR ) : mOp = MPI_BOR ; break;
        case(OP_BAND) : mOp = MPI_BAND; break;
        default      : check(-1, DMESG("No such MPI operation defined"));
    }
    return mOp;
}

#endif

#ifdef GKC_PARALLEL_MPI
MPI_Datatype Parallel::getMPIDataType(const std::type_info &T) {
    MPI_Datatype type=0;
    if     (T == typeid(double)) type = MPI_DOUBLE;
    else if(T == typeid(int   )) type = MPI_INT;
    else if(T == typeid(cmplxd)) type = MPI_DOUBLE_COMPLEX;
//    else if(T == typeid(bool  )) type = MPI_BOOL;
    else check(-1, DMESG("Such type is not defined"));
    
   return type;
}
#endif




Array2i Parallel::getProcessDomain(int rank) {
   int domain[6][3] = { { NxLD, NxLlD, NxLuD },
                        { NkyLD, NkyLlD, NkyLuD },
                        { NzLD, NzLlD, NzLuD },
                        { NvLD, NvLlD, NvLuD },
                        { NmLD, NmLlD, NmLuD },
                        { NsLD, NsLlD, NsLuD } };
    Array2i domain_out(Range(0,5), Range(0,2));

       // cannot use enum type here ... why !! DIR_X,...
#ifdef GKC_PARALLEL_MPI
                    MPI_Sendrecv(domain,  18, MPI_INT, master_rank,   99,   domain_out.data(),   18, MPI_INT, rank, 99, Comm[DIR_ALL], MPI_STATUS_IGNORE);
#endif
  return domain_out;
     }





Array1i Parallel::getAutoDecomposition(int numCPU) {
    Array1i A;
    A.resize(Range(0,5));
    A = 1;
    return A;

};
  

bool Parallel::checkValidDecomposition(Setup *setup, Array1i decomposition) {


   // Check basic decomposition sizes
   if( decomposition(DIR_X) > setup->get("Grid.Nx", 1)) check(-1, DMESG("Decomposition in x bigger than Nx"));
// no need to check y-decomposition (OpenMP parallelization)   
   if( decomposition(DIR_Y) > setup->get("Grid.Ny", 1)) check(-1, DMESG("Decomposition in y bigger than Ny"));
   if( decomposition(DIR_Z) > setup->get("Grid.Nz", 1)) check(-1, DMESG("Decomposition in z bigger than Nz"));
   if( decomposition(DIR_V) > setup->get("Grid.Nv", 1)) check(-1, DMESG("Decomposition in v bigger than Nv"));
   if( decomposition(DIR_M) > setup->get("Grid.Nm", 1)) check(-1, DMESG("Decomposition in m bigger than Nm"));
   if( decomposition(DIR_S) > setup->get("Grid.Ns", 1)) check(-1, DMESG("Decomposition in s bigger than Ns"));
   
   // Simple Check if reasonable values are provided for decomposition (only MPI proceeses)
   const int pNs = setup->get("Grid.Ns", 1 );
   //if((product(decomposition) != numThreads * numProcesses) && (myRank == 0)) check(-1, DMESG("Decomposition and number of processors are not equal"));
   if((product(decomposition) != numProcesses) && (myRank == 0)) check(-1, DMESG("Decomposition and number of processors are not equal"));
   if(((pNs %   decomposition(DIR_S)) != 0    ) && (myRank == 0)) check(-1, DMESG("Decomposition in s have to be modulo of the total number"));


    return GKC_SUCCESS;
};


void Parallel::print(std::string message) {

    if(myRank == 0) std::cout << message;

}


void Parallel::printOn(ostream &output) const {
       output << "Parallel   | ";
		
      if      (useOpenMP) output << " OpenMP (Threads)   : " << Num2String(numThreads)  << std::endl;
      if      (useOpenCL) output << " OpenCL (Threads)   : " << Num2String(numGPUs   )  << std::endl;
      if      (useMPI   ) {
	    output <<  " MPI (processes) : " << Num2String(numProcesses)  << " Decompostion : ";
            if (decomposition(DIR_X) == 0) output <<  "Automatic" << std::endl;
	    else
                     output << "X(" << decomposition(DIR_X) << ")  Y(" << decomposition(DIR_Y) << ") "
                            << "Z(" << decomposition(DIR_Z) << ")  V(" << decomposition(DIR_V) << ") "
                            << "M(" << decomposition(DIR_M) << ")  S(" << decomposition(DIR_S) << ") " << std::endl;
   
       } else output << std::endl;
     

#ifdef GKC_PARALLEL_MPI
/* 
//if(setup->flags & GKC_VERBOSE) {
                for(int rank = 0; rank < parallel.numProcesses; rank++) {
                    Array2i domain = parallel.getProcessDomain(rank);
                    if(parallel.myRank == 0) infoStream << 
                      
                      "Rank (" << rank<<")  | " <<
                        "X(" << domain(0, 0) << " " << domain(0, 1) << "-" << domain(0, 2) << ") " <<  
                        "Y(" << domain(1, 0) << " " << domain(1, 1) << "-" << domain(1, 2) << ") " <<  
                        "Z(" << domain(2, 0) << " " << domain(2, 1) << "-" << domain(2, 2) << ") " <<  
                        "V(" << domain(3, 0) << " " << domain(3, 1) << "-" << domain(3, 2) << ") " <<  
                        "M(" << domain(4, 0) << " " << domain(4, 1) << "-" << domain(4, 2) << ") " <<  
                        "S(" << domain(5, 0) << " " << domain(5, 1) << "-" << domain(5, 2) << ") " <<  
                       std::endl; 
                      
//                }
}
 * */
#endif // GKC_PARALLEL_MPI
    };
