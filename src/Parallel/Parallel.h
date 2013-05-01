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


#include <mpi.h>

// use enum class !
enum class Op {null = 0, sum=1, max=2, min=3, bor=4, band=5, land=6, lor=7};


extern int threadID;

/**
*  @brief Parallel Communication Interface
*
*  Current implementation supports using MPI-2
*  http://en.wikipedia.org/wiki/Message_Passing_Interface 
*  www.open-mpi.org
*
*  and OpenMP (3.1) http://www.openmp.org.
*
*  @todo move implementations to separate class
*
**/
class Parallel : public IfaceGKC {

/**
*  @brief  Class to simplify MPI handling, for each direction, such a class
*  
*  rank_l corresponds to lower rank in the dimension grid (..., myPosition - 1 , ... )
*  rank_u corresponds to upper rank in the dimension grid (..., myPosition + 1 , ... )
*
**/
struct NeighbourDir {
   MPI_Status  msg_status[4];  ///< Message status for 
   MPI_Request psf_msg_req[4]; ///< Message request id for phase space function
   int         rank_l,         ///< Lower rank
               rank_u;         ///< Upper rank
   int         psf_msg_tag[2], ///< Phase space message tag
               phi_msg_tag[2]; ///< Fields message tag
   NeighbourDir() {};
};

 public:
 
  /**
  *  @brief Decomposition of the code as 
  *         \f$ (x,y_k,z,v_\parallel, \mu, \sigma) \f$.
  *
  *  Note that $y_k$ decomposition is not supported
  *  by MPI but by OpenMP implementation
  *
  **/
  int decomposition[6];
   
  /**
  *  @brief Decomposition coordinate of process as 
  *         \f$ (x,y_k,z,v_\parallel, \mu, \sigma) \f$.
  *
  *  Coord[dir] for dir > DIR_Z coordinate is defined as e.g.
  *  Coord[DIR_VM] = (Coord[DIR_V] == 0) && (Coord[DIR_M] == 0).
  *  Warning : dis is not defined for all values of dir
  *
  **/
  int Coord[DIR_SIZE];
   
  int myRank;        ///< rank-id of process
  //    threadID;      ///< Number/Id of (OpenMP) thread  / possible to use extern thread_local ?! 

  int master_rank;   ///< rank-id of master process (usually 0)
   
  int numThreads,    ///< total number of OpenMP threads
      numProcesses;  ///< total number of MPI processes
   
  bool useOpenMP,    ///< true if compiled with OpenMP support
       useMPI;       ///< true if compiled with MPI support

  NeighbourDir Talk[DIR_S+1]; ///< Communication struct for Vlasov & Fields boundary

  MPI_Comm Comm[DIR_SIZE]; ///< MPI Communicators for Dir directions (warning not all implemented)
  int dirMaster[DIR_SIZE]; ///< Master process in direction dir (usually the one with Coord[dir] == 0)
  
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
  *   @todo should be possible to make constexpr
  *
  **/
  MPI_Datatype getMPIDataType(const std::type_info &T);

  Parallel(Setup *setup);
  virtual ~Parallel();
   
  /**
  *    @brief updates field boundaries for X,Z direction
  * 
  *    @todo how to make this allowed to be called only from Fields ?!
  *
  *
  *    @params  SendXl  Boundary values to send to Coord(dir +1)
  *    @params  SendXu  Boundary values to send to Coord(dir +1)
  *    @params  RecvXl  Boundary values to send to Coord(dir +1)
  *    @params  RecvXu  Boundary values to send to Coord(dir +1)
  *    @params  num_X   Number of elements for X-boundary
  *    @params  SendZl  Boundary values to send to Coord(dir +1)
  *    @params  SendZu  Boundary values to send to Coord(dir +1)
  *    @params  RecvZl  Boundary values to send to Coord(dir +1)
  *    @params  RecvZu  Boundary values to send to Coord(dir +1)
  *    @params  num_Z   Number of elements for Z-boundary
  *
  **/
  void updateBoundaryFields(CComplex *SendXl, CComplex *SendXu, CComplex *RecvXl, CComplex *RecvXu, int num_X,
                            CComplex *SendZl, CComplex *SendZu, CComplex *RecvZl, CComplex *RecvZu, int num_Z); 

   
   
  /**
  *    @brief updates Vlasov boundaries for direction dir
  * 
  *    @todo how to make this allowed to be called only from Vlasov ?!
  *
  *    This function uses non-blocking boundary conditions. As
  *    to calculate the field values a boundary points of the f
  *    are not required. However, once boundary values a required,
  *    we have to call updateBoundaryBarrier to make sure that
  *    boundary operations completed.
  *
  *    @params  Sendu  Boundary values to send to Coord(dir +1)
  *    @params  Sendl  Boundary values to send to Coord(dir +1)
  *    @params  Recvu  Boundary values to receive from Coord(dir +1)
  *    @params  Recvl  Boundary values to receive from Coord(dir +1)
  *    @params  num    Number of elements for dir-boundary
  *    @params  dir    Direction to send data
  *
  **/
  void updateBoundaryVlasov(CComplex *Sendu, CComplex *Sendl, CComplex *Recvu, CComplex  *Recvl, int num, int dir);

  /**
  *   
  *   @brief Barrier to assure that communication barriers completed.
  *
  *   @todo this is not intuitive, (use e.g. dir as input parameters)  
  *
  **/ 
  void updateBoundaryVlasovBarrier();
  
  /**
  *
  *    @brief bcast value to other process in direction dir
  *
  *    @param *A  pointer to array to send
  *    @param dir direction to send data
  *    @param num Number of data points
  *
  **/
  template<typename T> void bcast(T *A, int dir, int num) 
  {
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
  **/
  template<class T> T bcast(const T x, bool isRoot, int dir=DIR_ALL) 
  {

    // Notify all process who is root (is there a simpler way ?), take care it fails 
    // horribly if there is more than one root, (note 0 is master process also valid)
    // int master_rank = reduce(isRoot ? myRank : 0, Op::sum, dir);
    master_rank = 0;

#ifdef GKC_PARALLEL_MPI
    if(dir <= DIR_S) if(decomposition[dir] == 1) return x;
    // rank differ between communicators (?)
    check(MPI_Bcast((void *) &x, 1, getMPIDataType(typeid(T)), master_rank, Comm[dir]), DMESG("MPI_Bcast")); 
#endif
    return x; 
  }
   

  /**
  *   @brief Allreduce over direction dir
  *
  **/
  template<class T>  T reduce(T x, Op op, int dir=DIR_ALL, bool allreduce=true)
  {
#ifdef GKC_PARALLEL_MPI
    T global_dValue;

    check(MPI_Allreduce(&x, &global_dValue, 1, getMPIDataType(typeid(T)), getMPIOp(op), Comm[dir]), DMESG("MPI_Reduce")); 
    return global_dValue; 
#endif
    return x; 
  }
  
  /**
  *   @brief Allreduce over direction dir
  *
  *
  *
  **/
  template<class T>  void  reduce(T *A, Op op, int dir, int num, bool allreduce=true)
  {
    // No need to reduce if not decomposed in corresponding direction 
    // (should be handled by MPI library if in-place or ?
    if( (dir <= DIR_S) && (decomposition[dir] == 1)) return;

    // make range check for debugging
    if( (dir >= DIR_SIZE) ) check(-1, DMESG("Reduce : Invalid direction"));
  
    auto mpi_type = getMPIDataType(typeid(T)); 

#ifdef GKC_PARALLEL_MPI
    if(allreduce == true)
        check(MPI_Allreduce(MPI_IN_PLACE, A, num, mpi_type, getMPIOp(op), Comm[dir]), DMESG("MPI_AllReduce")); 
    else 
        check(MPI_Reduce((Coord[dir] == 0) ? MPI_IN_PLACE : A, A, num, mpi_type, getMPIOp(op), dirMaster[dir], Comm[dir]), DMESG("MPI_Reduce"   )); 
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
  *   @brief determines recommended decomposition 
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


  virtual void printProcessID(); 
  
  virtual void setThreadID(); 

  /**
  *  @brief initializes data output
  *
  *  Writes decomposition information, as well as version number
  *  underlying libraries
  *
  **/
  virtual void initData(Setup *setup, FileIO *fileIO);
 protected:

  virtual void printOn(std::ostream &output) const ;
  virtual void writeData(Timing *timing) {};
  virtual void closeData() {};
};

#endif // __SERIAL_H
