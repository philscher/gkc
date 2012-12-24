/*
 * =====================================================================================
 *
 *       Filename: Tools.h
 *
 *    Description: Loose collection of STATIC functions.
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __TOOLS_H_
#define __TOOLS_H_

/**
*
*   @brief collection of various helper functions
*
*
*
*
**/
class Tools
{

 public:

  /** 
  *     Create an integrer squence of [ 0, 1, 2, 3, 4, 5, 6 ...., N-1 ].
  *
  *     This squence is split around the number of CPUs in the communicator.
  *     
  *     e.g. for 3 CPUs and stop = 10, we will get a distribution like
  *
  *    CPU0 = [ 0, 1, 2, 3]
  *    CPU1 = [ 4, 5, 6]
  *    CPU2 = [ 7, 8, 9]
  *
  *    right now only MPI workers ( processes) are supported.
  *
  *
  *    ToDo : Is A more template based approach which also supports double etc useful ?
  *         
  *    Note : Don't expcet this to be performant (du eot creation of std::vector object)
  *
  *
  *
  **/
  static std::vector<int> ParallelRange(const int stop, Parallel *parallel, int dir) 
  {

    int Nworkers = parallel->getNumberOfWorkers(dir);
    int w_id     = parallel->getWorkerID(dir) ;
       
    // create list   
    int start = 0;

    int len   = stop - start;
    int chunk = len / Nworkers;
    
    int r = len % Nworkers;

    int w_start, w_stop;

    // what is r ? needs more documentation
    if (w_id < r) {
         w_start = start + (chunk + 1) * w_id;
         w_stop   = w_start + chunk;
    } else {
         w_start =   start + (chunk + 1) * r + chunk * (w_id - r);
         w_stop = w_start + chunk - 1;
    }

    std::vector<int> loop;// = { 1, 2, 3, 4, 5, 6, 7 };

    // create list from worker
    for(int n = w_start; n <= w_stop; n++) loop.push_back(n);

    return loop;
  };

};

#endif // __GKC_TOOLS_H__
