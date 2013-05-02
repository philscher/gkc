/*
 * =====================================================================================
 *
 *       Filename:  allocate.h
 *
 *    Description:  Allocates memory according to requested 
 *                  multi-dimensional range. Supports memory-align
 *                  allocation.
 *
 *                  Part of nct (numerical computing toolkit)
 *                  Contribute at https://github.com/philscher/nct
 *
 *         Author:  Paul Hilscher (c) 2012
 *
 *       License :  BSD
 *
 * =====================================================================================
 */

#ifndef __NCT_ARRAY_
#define __NCT_ARRAY_


#include <allocate.h>


//namespace nct { /// use better nct for numerical computing toolkit


#include<typeinfo>

/**
*    @brief Allocate four dimensional array using provided ranges
*
**/
template <class T> class Array4
{
  nct::allocate alloc; 
   
  nct::Range R0, R1, R2, R3;
  T *ptr_data0;
   
 public:  

  Array4() { };
   
  Array4(nct::Range R0, nct::Range R1, nct::Range R2, nct::Range R3) 
  { 
     // allocate array
   
     alloc = nct::allocate(R0, R1, R2, R3);
     alloc(&ptr_data0);
   
  };
   
 ~Array4()
  {
      // de-allocate
  };

  T& operator()(const int x0, const int x1, const int x2, const int x3)
  {

     // calculate index offset
     const int idx = 
             x0 * (R1.Num() * R2.Num() * R3.Num()) 
           + x1 * (R2.Num() * R3.Num()) 
           + x2 * (R3.Num()) 
           + x3 ;

       // check bounds
       
       return ptr_data0[idx];

  };
};
         
template <class T> class Array2
{
  nct::allocate alloc; 
   
  nct::Range R0, R1;
  T *ptr_data0;

 public:   
   
  Array2() { };
   
  Array2(nct::Range R0, nct::Range R1) 
  { 
  
     // allocate array
     alloc = nct::allocate(R0, R1);
     alloc(&ptr_data0);
   
  };
     
  T& operator()(const int x0, const int x1)
  {
     // calculate index offset
       
       const int idx = 
           + x0 * (R1.Num()) 
           + x1 ;

       // check bounds
       return ptr_data0[idx];
  };
   
  void operator=(T A)
  {
       T* p = alloc.data(&ptr_data0);

       for(int idx=0; idx<alloc.Num;idx++) p[idx] = T;
  };
};

template <class T> class Array1
{
  nct::allocate alloc; 
   
  nct::Range R0;
  T *ptr_data0;

 public:   
  
  Array1() { };
   
  Array1(nct::Range R0) 
  { 
     // allocate array
     alloc = nct::allocate(R0);
     alloc(&ptr_data0);
  };
     
  T& operator()(const int x0)
  {
     // calculate index offset
       const int idx = x0 ;

       return ptr_data0[idx];
  };
   
   void operator=(T A)
   {
       T* p = alloc.data(&ptr_data0);

       for(int idx=0; idx<alloc.Num;idx++) p[idx] = T;
   };

};

#endif //__NCT_ARRAY_
