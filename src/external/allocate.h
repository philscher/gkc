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

#ifndef __NCT_ALLOCATE_
#define __NCT_ALLOCATE_


#include <stdlib.h>
#include <iostream>
#include <vector>
#include <stack>

namespace nct { /// use better nct for numerical computing toolkit


/**
*   @brief simple range class
*   
*   @todo add copy constructor
*
**/
class Range
{

   int Start,  ///< Starting point
       Length; ///< Length of range


  public:
   Range(int _Start=0, int _Length=0) : Start(_Start), Length(_Length) {};

   /**
   *   @brief get number of elements inside the range
   *
   *   @return number of elements
   **/
   int Num() const { return Length; };

   /**
   *   @brief get starting point of range
   *
   *   @return starting point of range
   *
   **/
   int Off() const { return Start ; };

   void setRange(int _Start, int _Length) 
   {
      Start  = _Start;
      Length = _Length;
   };
};


enum alloc_flags  { 
               MA  =1      ,    ///< do use Memory Aligned (default)
               SET_ZERO=2  ,    ///< do not set to Zero 
               DEALLOC=4   ,    ///< do not deallocate 
               USE_DEFAULT = 7, ///< use default values 
             };


/**
*    @brief Allocating multi-dimensional arrays based on ranges
*
*    Currently up to 6-dimensional arrays are supported.
*
*    @note
*
*          Ranges offset can be positive. That means the pointer 
*          we have to artifically calculate the negative offset 
*          to compensate for this offset. 
*          However, it may happen that the negative offset is large
*          so that the final pointer address is negative.
*          
*          However, pointers on x86-64 are signed and symmetric around
*          0. Thus this offset calculation is valud. See
*
*          http://stackoverflow.com/questions/3304795/can-a-pointer-address-ever-be-negative
*
*     @rename R0, appropriatly first dimension is always continuous (access with stride = 1)
*
*
**/
#include<typeinfo>

class ArrayBase
{
   protected:
//     int Num;
//     int sizeOfType;
     
//     const std::type_info *typeInfo;
};


class allocate : public ArrayBase 
{
  int Num, ///< Number of Elements 
      Off; ///< Offset from p[0] to p[n] = first element

  int flags; // Deallocate all arrays after usage
  
  std::stack<void *> ptr_stack;

  public:

   int getNum() const { return Num; };
   int getOff() const { return Off; };

   allocate(int _flags = SET_ZERO | DEALLOC) : flags(_flags)
   {
      Num = 0;
      Off = 0;
   };


   /**
   *    @brief returns pointer with offset removed
   *
   *    Points to first element
   *
   **/
   template<class T> T* data(T *g0)
   {
         T *g = g0 + Off;
         return g; 
   };
   
   /**
   *
   *   @brief returns pointer which included offset
   *
   **/ 
   template<class T> T* zero(T *g)
   {
         T *g0 = g - Off;
         return g0; 
   };

   

   /**
   *    @brief calculated offset for one dimensional arrays.
   *
   *    @param R0         Range of first dimension
   *    @param user_flags Set additionally user flags
   *
   **/ 
   allocate(Range R0, alloc_flags user_flags = alloc_flags::USE_DEFAULT)             
   {
        flags &= user_flags;

        // get total array size
        Num = R0.Num() ; 
       
        // calculate offset to p[0]
        Off = R0.Off(); 
     
   };

   /**
   *    @brief calculated offset for two dimensional arrays.
   *
   *    @param R0         Range of first  dimension
   *    @param R1         Range of second dimension
   *
   *    @param user_flags Set additionally user flags
   *
   **/ 
   allocate(Range R0, Range R1, 
            alloc_flags user_flags = alloc_flags::USE_DEFAULT)             
   { 
        // get total array size
        Num = R0.Num() * R1.Num(); 
        
        // calculate offset to p[0][0]
        Off = 
              R0.Off() * (R1.Num()) 
            + R1.Off() ;
   };
   
   
   /**
   *    @brief calculated offset for three dimensional arrays.
   *
   *    Note : R0 has largest stride
   *           R2 is  continuous 
   *
   *    @param R0         Range of third  dimension
   *           R1         Range of second dimension
   *           R2         Range of first  dimension
   *
   *    @param user_flags Set additionally user flags
   *
   **/ 
   allocate(Range R0, Range R1, Range R2,
            alloc_flags user_flags = alloc_flags::USE_DEFAULT)             
   {
     
        // get total array size
        Num = R2.Num() * R1.Num() * R0.Num();

        // calculate offset to p[0][0][0]
        Off = 
          R0.Off() * (R1.Num() * R2.Num()) 
        + R1.Off() * (R2.Num()) 
        + R2.Off() ;

   };

   /**
   *    @brief calculated offset for four dimensional arrays.
   *
   *    @param R0 range in third  dimension
   *           R1 Range in second dimension
   *           R2 Range in first  dimension
   *           R3 Range in first  dimension
   *    @param user_flags Set additionally user flags
   *
   **/ 
   allocate(Range R0, Range R1, Range R2, Range R3,
            alloc_flags user_flags = alloc_flags::USE_DEFAULT)             
   {
       
       Num = R3.Num() * R2.Num() * R1.Num() * R0.Num();

       Off = 
          R0.Off() * (R1.Num() * R2.Num() * R3.Num()) 
        + R1.Off() * (R2.Num() * R3.Num()) 
        + R2.Off() * (R3.Num()) 
        + R3.Off() ;

   };
   
   /**
   *    @brief calculated offset for five dimensional arrays.
   *
   *    Note : R0 has largest stride
   *           R2 is  continuous 
   *
   *    @param R0 range in third  dimension
   *           R1 Range in second dimension
   *           R2 Range in first  dimension
   *           R2 Range in first  dimension
   *           R2 Range in first  dimension
   *           R2 Range in first  dimension
   *    @param user_flags Set additionally user flags
   *
   **/ 
   allocate(Range R0, Range R1, Range R2, Range R3,
            Range R4,
            alloc_flags user_flags = alloc_flags::USE_DEFAULT)             
   {
       
       Num = R4.Num() * R3.Num() * R2.Num() * 
             R1.Num() * R0.Num();
       
       Off = 
          R0.Off() * (R1.Num() * R2.Num() * R3.Num() * R4.Num()) 
        + R1.Off() * (R2.Num() * R3.Num() * R4.Num()) 
        + R2.Off() * (R3.Num() * R4.Num()) 
        + R3.Off() * (R4.Num()) 
        + R4.Off() ;

   };

   /**
   *    @brief calculated offset for six dimensional arrays.
   *
   *    Note : R0 has largest stride
   *           R2 is  continuous 
   *
   *    @param R0 range in third  dimension
   *           R1 Range in second dimension
   *           R2 Range in first  dimension
   *    @param user_flags Set additionally user flags
   *
   **/ 
   allocate(Range R0, Range R1, Range R2, Range R3,
            Range R4, Range R5,
            alloc_flags user_flags = alloc_flags::USE_DEFAULT)             
   {
       // get total array size
       Num = R5.Num() * R4.Num() * R3.Num() * 
             R2.Num() * R1.Num() * R0.Num();
       
       // calculate offset to p[0][0][0][0][0][0]
       Off = 
          R0.Off() * (R1.Num() * R2.Num() * R3.Num() * R4.Num() * R5.Num()) 
        + R1.Off() * (R2.Num() * R3.Num() * R4.Num() * R5.Num()) 
        + R2.Off() * (R3.Num() * R4.Num() * R5.Num()) 
        + R3.Off() * (R4.Num() * R5.Num()) 
        + R4.Off() * (R5.Num()) 
        + R5.Off() ;

   };

   
    /**
    *    @brief deallocate arrays
    *
    *    All arrays allocated with this class will get release
    *    once it's instance is destroyed if DEALLOC (default)
    *    flag is set.
    *
    **/ 
    ~allocate()
    {
         if(flags & alloc_flags::DEALLOC) { 
/* 
            // Deallocate all arrays
            while (!ptr_stack.empty()) {
           
              // release top element
              void *ptr = ptr_stack.top();
               
              //if() free(ptr)
              _mm_free(ptr);
              // remove element on top
              ptr_stack.pop();
            }
 * */
         }
    };

    /**
    *     @brief allocate memory space for function.
    *
    *     In order to increase performance, we should align the data. 
    *     The best alignment can be found by using ...
    *
    *
    **/
    template<class T> allocate& operator()(T **g)
    {
         
         //*g = ((T *) malloc(Num * sizeof(T)));
         *g = ((T *) _mm_malloc(Num * sizeof(T), 64));

           ptr_stack.push((void *) *g);
        
         // set to zero
         if(flags & SET_ZERO) {
            for(int n=0; n < Num; n++) (*g)[n] = T(0);
         }

         // Take care, pointer arithmetic is typed, only char* is 1 Byte !!!
         // Substract offset to calculate p[0][0]....
         *g = *g - Off;
    
         return *this;
    };
 
    /**
    *     @brief allocate memory space for function.
    *
    *     In order to increase performance, we should align the data. 
    *     The best alignment can be found by using ...
    *
    *
    **/
    template<class T> allocate& operator()(T **g0, T **g1)
    {
         operator()(g0);
         operator()(g1);

         return *this;
    };
 
    /**
    *     @brief allocate memory space for function.
    *
    *     In order to increase performance, we should align the data. 
    *     The best alignment can be found by using ...
    *
    *
    **/
    template<class T> allocate& operator()(T **g0, T **g1, T **g2)
    {
         operator()(g0);
         operator()(g1);
         operator()(g2);
         
         return *this;
    };
 
    /**
    *     @brief allocate memory space for function.
    *
    *     In order to increase performance, we should align the data. 
    *     The best alignment can be found by using ...
    *
    *
    **/
    template<class T> allocate& operator()(T **g0, T **g1, T **g2, T **g3)
    {
         operator()(g0);
         operator()(g1);
         operator()(g2);
         operator()(g3);
         
         return *this;
    }
 
    /**
    *     @brief allocate memory space for function.
    *
    *     In order to increase performance, we should align the data. 
    *     The best alignment can be found by using ...
    *
    *
    **/
    template<class T> allocate& operator()(T **g0, T **g1, T **g2, T **g3, T **g4)
    {
         operator()(g0);
         operator()(g1);
         operator()(g2);
         operator()(g3);
         operator()(g4);
         
         return *this;
    }
 
    /**
    *     @brief allocate memory space for function.
    *
    *     In order to increase performance, we should align the data. 
    *     The best alignment can be found by using ...
    *
    *
    **/
    template<class T> allocate& operator()(T **g0, T **g1, T **g2, T **g3, T **g4, T **g5)
    {
         operator()(g0);
         operator()(g1);
         operator()(g2);
         operator()(g3);
         operator()(g4);
         operator()(g5);
         
         return *this;
    }

    template<class T> allocate& operator()(T **g0, T **g1, T **g2, T **g3, T **g4, T **g5, T**g6)
    {
         operator()(g0);
         operator()(g1);
         operator()(g2);
         operator()(g3);
         operator()(g4);
         operator()(g5);
         operator()(g6);
         
         return *this;
    }

 };

} // namespace nct 


#endif // __NCT_ALLOCATE_
