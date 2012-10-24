/*
 * =====================================================================================
 *
 *       Filename: Parallel.h
 *
 *    Description: Parallelization functions declaration.
 *
 *         Author: Paul P. Hilscher (2010), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __PARALELL_H
#define __PARALELL_H


#include "Global.h"
#include "Setup.h"

#ifdef GKC_PARALLEL_OPENMP
#include <omp.h>
#endif

// use enum class !
enum class Op {OP_NULL = 0, SUM=1, MAX=2, MIN=3, BOR=4, BAND=5, LAND=6, LOR=7};

#include <mpi.h>


/**
*  @brief Interface for MPI and OpenMP
*  
*
**/
class Parallel : public IfaceGKC {

/**
*  @brief  Class to simplify MPI handling, for every direction we define such a class
*  
*  rank_l correspindes to lower rank in the dimension grid (..., myPosition - 1 , ... )
*  rank_u correspindes to upper rank in the dimension grid (..., myPosition + 1 , ... )
*
**/
struct NeighbourDir {
   MPI_Status  msg_status[4];
   MPI_Request psf_msg_request[4];
   int         rank_l,
               rank_u;
   int         psf_msg_tag[2], phi_msg_tag[2];
   NeighbourDir() {};
};

  public:
 

   /**
   *  @brief MPI Decompostion of the code (X,Y,Z,V,M,S)
   *
   **/
   int decomposition[6];
   
   /**
   *   @brief Coordinate of the process in (X,Y,Z,V,M,S)
   *
   **/
   int Coord[DIR_SIZE];
   
   //! hold the (world) rank of the process
   int myRank; // substitute myRank with mpi_rank
   int master_rank;
   
   //! the total number of threads 
   int numThreads, numProcesses;
   //! switches if OpenMP, MPI is during compule time
   bool useOpenMP, useMPI;

   // Some MPI specific stuff
   NeighbourDir Talk[DIR_S+1];
   MPI_Status stat; 

   MPI_Comm Comm[DIR_SIZE];
   int dirMaster[DIR_SIZE];
  
   /**
   *   @brief barrier
   *
   **/
   void barrier(int dir=DIR_ALL);


   /**
   *   @brief get corresponding MPI operation from data type
   *
   **/
   MPI_Op   getMPIOp(Op op);

   /**
   *   @brief get MPI data type from  C++ type
   *
   **/
   MPI_Datatype getMPIDataType(const std::type_info &T);

   Parallel(Setup *setup);
   virtual ~Parallel();
   
   /**
   *    @brief updates boundaries in X,Z direction
   *
   *
   **/
   void  updateBoundaryFields(CComplex *SendXl, CComplex *SendXu, CComplex *RecvXl, CComplex *RecvXu, int numElements_X,
                              CComplex *SendZl, CComplex *SendZu, CComplex *RecvZl, CComplex *RecvZu, int numElements_Z); 

   
   
   template<class T> void updateBoundaryVlasov(T *Sendu, T *Sendl, T *Recvu, T  *Recvl, int numElements, int dir) {
        
#ifdef GKC_PARALLEL_MPI
     MPI_Irecv(Recvl, numElements, getMPIDataType(typeid(T)), Talk[dir].rank_l, Talk[dir].psf_msg_tag[0], Comm[DIR_ALL], &Talk[dir].psf_msg_request[1]);
     MPI_Isend(Sendu, numElements, getMPIDataType(typeid(T)), Talk[dir].rank_u, Talk[dir].psf_msg_tag[0], Comm[DIR_ALL], &Talk[dir].psf_msg_request[0]);
     
     if(Talk[dir].rank_u == MPI_PROC_NULL)   Recvl[0:numElements] = 0.e0;
     
     MPI_Irecv(Recvu, numElements, getMPIDataType(typeid(T)), Talk[dir].rank_u, Talk[dir].psf_msg_tag[1], Comm[DIR_ALL], &Talk[dir].psf_msg_request[3]); 
     MPI_Isend(Sendl, numElements, getMPIDataType(typeid(T)), Talk[dir].rank_l, Talk[dir].psf_msg_tag[1], Comm[DIR_ALL], &Talk[dir].psf_msg_request[2]);
     if(Talk[dir].rank_l == MPI_PROC_NULL)   Recvu[0:numElements] = 0.e0;
#endif // GKC_PARALLEL_MPI
     return;

   };


   void updateNeighboursBarrier();
   
   template<typename T> void send(T *A, int dir, int num) {
#ifdef GKC_PARALLEL_MPI
     if(dir <= DIR_S) if(decomposition[dir] == 1) return;
     // rank differ between communicators, so we should take care of rankFFTMaster
     check(MPI_Bcast(A,num, getMPIDataType(typeid(T)), dirMaster[dir], Comm[dir]), DMESG("MPI_Bcast")); 
#endif
     return;
   };



   /**
   *   @brief sends scalar data to other CPU 
   *
   *   @todo rename to bcast
   *
   **/
   template<class T> void send(T &x, bool isRoot, int dir=DIR_ALL) 
   {

     // Notify all process who is root (is there a simpler way ?), take care it fails 
     // horribly if there is more than one root, (note 0 is master process also valid)
     int master_rank = collect(isRoot ? myRank : 0, Op::SUM, dir);
     master_rank = 0;

#ifdef GKC_PARALLEL_MPI
     if(dir <= DIR_S) if(decomposition[dir] == 1) return;
     // rank differ between communicators, so we should take care of rankFFTMaster
     check(MPI_Bcast(&x, 1, getMPIDataType(typeid(T)), master_rank, Comm[dir]), DMESG("MPI_Bcast")); 
#endif
     return; 
   }
   

   /**
   *   @brief Allreduce over direction dir
   *
   **/
   template<class T>  T  collect(T x, Op op, int dir=DIR_ALL, bool allreduce=true)
   {
#ifdef GKC_PARALLEL_MPI
     T global_dValue;

     // we need allreduce instead of reduce because H5TB need all process to have the same value
     MPI_Allreduce(&x, &global_dValue, 1, getMPIDataType(typeid(T)), getMPIOp(op), Comm[dir]);//, DMESG("MPI_Reduce"); 
     //check(MPI_Allreduce(&x, &global_dValue, N, getMPIDataType(typeid(T)), getMPIOp(op), Comm[dir]), DMESG("MPI_Reduce")); 
     return global_dValue; 
#endif
     return x; 
   }
   
   
   /**
   *   @brief Allreduce over direction dir
   *
   **/
   template<class T>  void  collect(T *A, Op op, int dir, int Num, bool allreduce=true)
   {
     //if(dir <= DIR_S) if(decomposition[dir] == 1) return A;
     
#ifdef GKC_PARALLEL_MPI
     // note, for MPI_Reduce only root process can specify MPI_IN_PLACE
     // we need allreduce instead of reduce because H5TB need all process to have the same value
        //MPI_Allreduce(MPI_IN_PLACE, A.data(), A.numElements(), getMPIDataType(typeid(T)), getMPIOp(op),                Comm[dir]), DMESG("MPI_Allreduce");
     if(allreduce == true)
        check(MPI_Allreduce(MPI_IN_PLACE, A, Num, getMPIDataType(typeid(T)), getMPIOp(op), Comm[dir]), DMESG("MPI_AllReduce")); 
     else 
        check(MPI_Reduce((myRank == dirMaster[dir]) ? MPI_IN_PLACE : A, A, Num, getMPIDataType(typeid(T)), getMPIOp(op), dirMaster[dir], Comm[dir]), DMESG("MPI_Reduce"   )); 
#endif
     return;
   }

   /**
   *   @brief gets total number of process in direction DIR
   *
   *   @param  dir see enum Direction
   *   @return     number of processes
   **/
   int getNumberOfWorkers(int dir);
    
   
   /**
   *   Please Document me
   **/
   int getWorkerID(int dir) ;


   /**
   *    @brief simple checks determine if decomposition is logical correct
   **/
   void   checkValidDecomposition(Setup *setup);

   /**
   *   @brief determines recommendend decomposition 
   *
   *   @param numCPU Total Number of available CPUs
   *   @return  Recommended decomposition
   *
   **/
   void getAutoDecomposition(int numCPU);
   
   /**
   *  @brief  Prints out string to terminal (only master process)
   *
   **/
   virtual void print(std::string message);

  protected:

   virtual void printOn(std::ostream &output) const ;
   virtual void initData(FileIO *fileIO) {};
   virtual void writeData(Timing *timing) {};
   virtual void closeData() {};
};

#endif // __SERIAL_H
