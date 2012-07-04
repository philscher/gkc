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

#define PARALLEL_OPENMP

#ifdef PARALLEL_OPENMP
#include <omp.h>
#endif

enum Operations {OP_NULL = 0, OP_SUM=1, OP_MAX=2, OP_MIN=3, OP_BOR=4, OP_BAND};

#ifdef GKC_PARALLEL_MPI

#include <mpi.h>

//! Class to simplify MPI handling, for every direction we define such a class
/*!
 *     rank_l correspindes to lower rank in the dimension grid (..., myPosition - 1 , ... )
 *     rank_u correspindes to upper rank in the dimension grid (..., myPosition + 1 , ... )
 *
 *
 */
struct NeighbourDir {
   MPI_Status  msg_status[4];
   MPI_Request psf_msg_request[4];
   int         rank_l,        rank_u;
   int         psf_msg_tag[2], phi_msg_tag[2];
   NeighbourDir() {};
};

#endif



//! Class handles parallelization with OpenMP, MPI and OpenCL (future)
/*!
 *  Parallel describes the communication between grid ...
 *  
 *  \image html mpi.png
 *
 */
class Parallel : public IfaceHelios {

  public:
    Array1i Coord;
    Array1i decomposition;
  //! hold the (world) rank of the process
  int myRank, master_process_id;
  int master_rank;
  bool isAutoDecomposed;
  //! the total number of threads 
  int numThreads, numProcesses, numGPUs;
  //! switches if OpenMP, MPI or OpenCL is during compule time
  bool useOpenMP, useMPI, useOpenCL;

#ifdef GKC_PARALLEL_MPI
  Array<NeighbourDir, 1> Talk;
   MPI_Status stat; 
//   int rankMyY0, rankMyX0, rankFFT1DMaster;

   MPI_Comm Comm[DIR_SIZE];
   int dirMaster[DIR_SIZE];
   
   void barrier(int dir=DIR_ALL) {
        MPI_Barrier(Comm[dir]);
   };

//   MPI_Comm getComm(int dir);
   MPI_Op   getMPIOp(int op);
MPI_Datatype getMPIDataType(const std::type_info &T);

#endif //Parallel_MPI
   Parallel(Setup *setup);
  virtual ~Parallel();
   
  
   int  updateNeighbours(Array6z  SendXl, Array6z  SendXu, Array6z  SendYl, Array6z  SendYu, Array6z SendZl, Array6z SendZu, 
                         Array6z  RecvXl, Array6z  RecvXu, Array6z  RecvYl, Array6z  RecvYu, Array6z RecvZl, Array6z RecvZu); 

   
   //* Sets Recv to 0.e0 if it should recv from a program with MPI_PROC_NULL (at the end of the domain with
   // no periodic boundary conditions. e.g. for velocity.
   template<typename T, int W> int updateNeighbours(Array<T,W>  Sendu,  Array<T,W>  Sendl,  Array<T,W>  Recvu, Array<T,W>  Recvl, int dir) {
        
#ifdef GKC_PARALLEL_MPI
        MPI_Irecv(Recvl.data(), Recvl.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_l, Talk(dir).psf_msg_tag[0], Comm[DIR_ALL], &Talk(dir).psf_msg_request[1]);
        MPI_Isend(Sendu.data(), Sendu.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_u, Talk(dir).psf_msg_tag[0], Comm[DIR_ALL], &Talk(dir).psf_msg_request[0]);
        if(Talk(dir).rank_u == MPI_PROC_NULL)   Recvl = 0.e0;
     
        MPI_Irecv(Recvu.data(), Recvu.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_u, Talk(dir).psf_msg_tag[1], Comm[DIR_ALL], &Talk(dir).psf_msg_request[3]); 
        MPI_Isend(Sendl.data(), Sendl.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_l, Talk(dir).psf_msg_tag[1], Comm[DIR_ALL], &Talk(dir).psf_msg_request[2]);
        if(Talk(dir).rank_l == MPI_PROC_NULL)   Recvu = 0.e0;
#endif // GKC_PARALLEL_MPI
      return HELIOS_SUCCESS;

   };
   int updateNeighboursBarrier();
   
   template<typename T, int W> int updateNeighbours(Array<T,W>  Sendu,  Array<T,W>  Sendl,  Array<T,W>  Recvu, Array<T,W>  Recvl, int dir, bool nonBlocking) {
#ifdef GKC_PARALLEL_MPI
        int msg_tag[] = { 4998 , 4999};
        MPI_Status  msg_status[2];

        MPI_Sendrecv(Sendu.data(), Sendu.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_u, msg_tag[1], 
                     //Recvl.data(), Recvl.umElements(), getMPIDataType(typeid(T)), Talk(dir).rank_l, msg_tag[1], Comm[dir],  &msg_status[0]);
                     Recvl.data(), Recvl.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_l, msg_tag[1], Comm[DIR_ALL],  &msg_status[0]);
        MPI_Sendrecv(Sendl.data(), Sendl.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_l, msg_tag[0], 
                     //Recvu.data(), Recvu.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_u, msg_tag[0], Comm[dir], &msg_status[1]);
                     Recvu.data(), Recvu.numElements(), getMPIDataType(typeid(T)), Talk(dir).rank_u, msg_tag[0], Comm[DIR_ALL], &msg_status[1]);
 
#endif // GKC_PARALLEL_MPI
      return HELIOS_SUCCESS;

   };
   
 
template<typename T, int W> int send(Array<T,W> A, int dir=DIR_ALL) {
#ifdef GKC_PARALLEL_MPI
     if(dir <= DIR_S) if(decomposition(dir) == 1) return HELIOS_SUCCESS;
    // rank differ between communicators, so we should take care of rankFFTMaster
     check(MPI_Bcast(A.data(), A.numElements(), getMPIDataType(typeid(T)), dirMaster[dir], Comm[dir]), DMESG("MPI_Bcast")); 
#endif
       return HELIOS_SUCCESS; 
}

template<class T> int send(T &x, bool isRoot, int dir=DIR_ALL) {

    // Notify all process who is root (is there a simpler way ?), take care it fails 
    // horribly if there is more than one root, (note 0 is master process also valid)
    int master_rank = collect(isRoot ? myRank : 0, OP_SUM, dir);
    master_rank = 0;

#ifdef GKC_PARALLEL_MPI
     if(dir <= DIR_S) if(decomposition(dir) == 1) return HELIOS_SUCCESS;
    // rank differ between communicators, so we should take care of rankFFTMaster
     check(MPI_Bcast(&x, 1, getMPIDataType(typeid(T)), master_rank, Comm[dir]), DMESG("MPI_Bcast")); 
#endif
       return HELIOS_SUCCESS; 
}
   
// Fix this
template<class T>  T  collect2(T x, int numElements=1, int op = OP_SUM, int dir=DIR_ALL, bool allreduce=true)
{
#ifdef GKC_PARALLEL_MPI
        T global_dValue;
       // we need allreduce instead of reduce because H5TB need all process to have the same value
       check(MPI_Allreduce(&x, &global_dValue, numElements, getMPIDataType(typeid(T)), getMPIOp(op), Comm[dir]), DMESG("MPI_Reduce")); 
       return global_dValue; 
#endif
       return x; 
}

template<class T>  T  collect(T x, int op = OP_SUM, int dir=DIR_ALL, bool allreduce=true)
{
#ifdef GKC_PARALLEL_MPI
        T global_dValue;
       // we need allreduce instead of reduce because H5TB need all process to have the same value
       check(MPI_Allreduce(&x, &global_dValue, 1, getMPIDataType(typeid(T)), getMPIOp(op), Comm[dir]), DMESG("MPI_Reduce")); 
       return global_dValue; 
#endif
       return x; 
}

template<typename T, int W> Array<T,W> collect(Array<T,W> A, int op=OP_SUM, int dir=DIR_ALL, bool allreduce=true) {
#ifdef GKC_PARALLEL_MPI
     // Return immediately if we don't decompose in this direction
     if(dir <= DIR_S) if(decomposition(dir) == 1) return A;
    
     if(allreduce == true)
        MPI_Allreduce(MPI_IN_PLACE, A.data(), A.numElements(), getMPIDataType(typeid(T)), getMPIOp(op),                Comm[dir]), DMESG("MPI_Allreduce");
 //       check(MPI_Allreduce(MPI_IN_PLACE, A.data(), A.numElements(), getMPIDataType(typeid(T)), getMPIOp(op),                Comm[dir]), DMESG("MPI_Allreduce")); 
     else    // note, for MPI_Reduce only root process can specify MPI_IN_PLACE
        check(MPI_Reduce   ((myRank == dirMaster[dir]) ? MPI_IN_PLACE : A.data(), A.data(), A.numElements(), getMPIDataType(typeid(T)), getMPIOp(op), dirMaster[dir], Comm[dir]), DMESG("MPI_Reduce"   )); 
#endif
     return A;
}

     Array2i getProcessDomain(int rank);

     /**
      *   Prints out string to terminal (only master process)
      *    
      *
      *
      * */
     virtual void print(std::string message);
     virtual void print(const char *message) { print(std::string(message));};
     /**
      *   Prints out stringstream to terminal (only master process)
      *    
      *
      *
      * */
     virtual void print(std::stringstream &message) { 
       std::string msg_str = message.str();
       print(msg_str); 
     }; 


    int getNumberOfWorkers(int dir) {
        int numWorkers = 0;
        MPI_Comm_size(Comm[dir],&numWorkers);
        return numWorkers;
    };
    
    int getWorkerID(int dir) {
        int rankWorker = 0;
        MPI_Comm_rank(Comm[dir],&rankWorker);
        return rankWorker;
    };

  // friend ostream& operator<<(ostream& output, const Parallel& parallel) {
   virtual void printOn(ostream &output) const {
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
//if(setup->flags & HELIOS_VERBOSE) {
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
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

bool    checkValidDecomposition(Setup *setup, Array1i decomposition);
Array1i getAutoDecomposition(int numCPU);

};

#endif // __SERIAL_H
