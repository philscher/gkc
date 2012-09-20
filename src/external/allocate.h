/*
 * =====================================================================================
 *
 *       Filename:  array_sub.cpp
 *
 *    Description:  i
 *
 *         Author:  Paul Hilscher (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <iostream>

namespace hpc {

class Range
{
   int Start, Length;

  public:
   Range(int _Start, int _Length) : Start(_Start), Length(_Length) {};

   int getNum() const { return Length; };
   int getOff() const { return Start ; };
   
   Range getOff(Range Rn)  { return Range(0, Rn.getOff() * Length + Start); };
   //Range getOff(Range Rn)  { return Range(0, Start * Rn.getNum() + Rn.getOff()); };

};
//http://stackoverflow.com/questions/3304795/can-a-pointer-address-ever-be-negative
// on x86-64 pointers are sign extended and symmetric around 0 , thus the
// offset calculation is valid.

class allocate 
{
  int N, Off;

  public:

   // calculate offsets
   allocate(Range R0)             { N = R0.getNum() ; 
                                     Off = R0.getOff(); 
                                   };
   allocate(Range R0, Range R1)  { N = R0.getNum() * R1.getNum(); 
                                     Off = R0.getOff(R1).getNum();
                                     std::cout << "Using 2D !!!!!!!!!!!!!!!" << std::endl;
                                   };
   allocate(Range R0, Range R1, Range R2) 
   {
     std::cout
       << "R0 : " << R0.getNum() << " - " << R0.getOff() << std::endl
       << "R1 : " << R1.getNum() << " - " << R1.getOff() << std::endl
       << "R2 : " << R2.getNum() << " - " << R2.getOff() << std::endl;
     N = R2.getNum() * R1.getNum() * R0.getNum();

     Off = 
       R0.getOff() * (R1.getNum() * R2.getNum()) 
     + R1.getOff() * (R2.getNum()) 
     + R2.getOff() ;
     std::cout << "Offset : " << Off  << std::endl;
     //Off = R0.getOff(R2).getOff(R1).getNum(); 
   };
   allocate(Range R0, Range R1, Range R2, Range R3) 
   {
     N = R0.getNum() * R1.getNum() * R2.getNum() * R3.getNum(); 
     Off = R0.getOff(R1).getOff(R2).getOff(R3).getOff(); 
   };
   
   allocate(Range R0, Range R1, Range R2, Range R3,
            Range R4) 
   {
     N = R0.getNum() * R1.getNum() * R2.getNum() * R3.getNum() *
         R4.getNum(); 
     Off = R0.getOff(R1).getOff(R2).getOff(R3).getOff(R4).getOff(); 
   };
   
   allocate(Range R0, Range R1, Range R2, Range R3,
            Range R4, Range R5) 
   {
     N = R0.getNum() * R1.getNum() * R2.getNum() * R3.getNum() *
         R4.getNum() * R5.getNum(); 
     Off = R0.getOff(R1).getOff(R2).getOff(R3).getOff(R4).getOff(R5).getOff(); 
   };


 // allocate variables

 template<class T> void operator()(T **g)
 {
//   abort();
//  using mm_malloc crashes the code, why ?! !
   // Take care, pointer arithmetic is typed, only char* is 1 Byte !!!
   //*g = (T *) (((char *) malloc(N * sizeof(T)))   - Off * sizeof(T));
   *g = ((T *) malloc(N * sizeof(T))) - Off;
 // Not working ...   *g = (T *) (((char *) _mm_malloc(N * sizeof(T),64))   - Off * sizeof(T));
 };
 
 template<class T> void operator()(T **g0, T **g1)
 {
   operator()(g0);
   operator()(g1);
 };
 
 template<class T> void operator()(T **g0, T **g1, T **g2)
 {
   operator()(g0);
   operator()(g1);
   operator()(g2);
};
 
 template<class T> void operator()(T **g0, T **g1, T **g2, T **g3)
 {
   operator()(g0);
   operator()(g1);
   operator()(g2);
   operator()(g3);
 }
 
 template<class T> void operator()(T **g0, T **g1, T **g2, T **g3, T **g4)
 {
   operator()(g0);
   operator()(g1);
   operator()(g2);
   operator()(g3);
   operator()(g4);
 }
 
 template<class T> void operator()(T **g0, T **g1, T **g2, T **g3, T **g4, T **g5)
 {
   operator()(g0);
   operator()(g1);
   operator()(g2);
   operator()(g3);
   operator()(g4);
   operator()(g5);
 }
};

//allocate(Range(NxLlD,Nx),Range(NyLlD, Ny))(&g);
} 
