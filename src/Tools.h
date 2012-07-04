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

class Tools
{
  public:
/* 
     std::list Range(const int stop) {

       std::list<int> L;
       if (stop > 0 ) for(int n = 0 ; n < stop; n++)  L.push_back(n);
     
     }

     }  
     std::list Range(const int start, const int stop, const int step) {
       std::list<int> range;

       if(start => stop

       the_list.push_back(1);
       the_list.push_back(-15);
       the_list.push_back(3);

            // create full list first
     }
 * */



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
      * */
     static std::vector<int> ParallelRange(const int stop, Parallel *parallel, int dir) {

        //	    int Nworkers = parallel->getNumberWorkers();
        int Nworkers = parallel->getNumberOfWorkers(dir);
        int w_id = parallel->getWorkerID(dir) ;
	    // create list	

        int start = 0;


        int len   = stop - start;
        int chunk = len / Nworkers;

        int r = len % Nworkers;
   
        int w_start, w_stop;


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

   //   for(auto n : loop) std::cout << n << " " ;
   //   std::cout << std::endl;

      return loop;
    };

};

#endif // __TOOLS_H_
