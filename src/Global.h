/*
 * =====================================================================================
 *
 *       Filename: Global.h
 *
 *    Description: Global variables/functions definitions
 *
 *         Author: Paul P. Hilscher (2009), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GLOBAL_H
#define __GLOBAL_H


#include <cilk/cilk.h>


#include<string>
#include<complex>
#include<iostream>
#include<typeinfo>
#include<iomanip>

#include "external/allocate.h"
#include "config.h"

///////////////////// Define data types ////////////////
typedef std::complex<double> Complex;  

typedef _Complex double CComplex;  
typedef double               Real   ;
#define _imag ((CComplex) (0.+1.j)) 
   
// align to cache-lines
typedef __declspec(align(64)) double     doubleAA;
typedef __declspec(align(64)) CComplex CComplexAA;

////////////////////////////////////////////////////////

extern "C" double cabs  (CComplex z);
extern "C" double creal (CComplex z);
extern "C" double cimag (CComplex z);
extern "C" double carg  (CComplex z);
extern "C" CComplex conj(CComplex z);
extern "C" CComplex cexp (CComplex z);




template<class T> __attribute__((vector)) inline T pow2(T x) { return x*x; };
template<class T> inline T pow3(T x) { return x*x*x; };
template<class T> inline T pow4(T x) { const T x2 =  (x*x); return x2*x2; };
template<class T> inline T pow5(T x) { const T x2 = x * x; return x2*x2*x; };
template<class T> inline T pow6(T x) { const T x2 = x * x; return x2*x2*x2; };
template<class T> inline T pow8(T x) { const T x2 = x * x; x2 *= x2; return x2*x2; };
//__declspec(vector) inline CComplex square(CComplex x) { return x*x; };
inline CComplex square(CComplex x) { return x*x; };

// for vectorization support, see 
#define simd_for _Pragma("simd") for


// DIR_SIZE has to be the number of elements, excluding DIR_SIZE
enum Dir : int {DIR_X=0, DIR_Y, DIR_Z, DIR_V, DIR_M, DIR_S, DIR_ALL, DIR_XYZ, DIR_VMS, DIR_MS, DIR_VM, DIR_XY, DIR_XYZVM, DIR_SIZE};
     
inline int check( int status, std::string file, int line, std::string error_text, bool doAbort=false)
{
  if(status == -1 ) {
            
    std::stringstream ss;
    
    ss << std::endl; 
    ss << "\033[0;m";
    ss << "\033[1;m"  <<"!.... " << file << "(" << line << ") " << error_text; 
    ss << "\033[0m" << std::endl << std::flush;

    std::cout << ss.str();
              
    // exit through abort so we can get stack trace
    if(doAbort==true) abort();
    else              exit(0);
    
  }
  
  return status;

}


// use assert instead !
#define DMESG(mesg)  std::string(__FILE__), __LINE__, std::string(mesg)
#define _D(mesg, value) std::cout << mesg << " " << value << std::endl << std::flush



/*************************  Global Variables   ********************************/

// DO NOT CHANGE THESE VARIABLES BY ANY MEANS OUTSIDE GRID
// OR CONFIG !!

// Gx[LlD], Gx[LuD], Gx[LuB], Gx[LuD] // should improve speed due to caching
extern int Nx, Nky, Nz, Nv, Nm, Ns;

extern int NxLlD, NxLuD, NxLlB, NxLuB; 
extern int NyLlD, NyLuD, NyLlB, NyLuB; 
extern int NzLlD, NzLuD, NzLlB, NzLuB; 
extern int NvLlD, NvLuD, NvLlB, NvLuB; 
extern int NmLlD, NmLuD, NmLlB, NmLuB; 
extern int NsLlD, NsLuD, NsLlB, NsLuB; 

extern int NxGlD, NxGuD, NxGlB, NxGuB; 
extern int NyGlD, NyGuD, NyGlB, NyGuB; 
extern int NzGlD, NzGuD, NzGlB, NzGuB; 
extern int NvGlD, NvGuD, NvGlB, NvGuB; 
extern int NmGlD, NmGuD, NmGlB, NmGuB; 
extern int NsGlD, NsGuD, NsGlB, NsGuB; 


extern int NkyGlD, NkyGuD, NkyGlB, NkyGuB; 
extern int NkyLlD, NkyLuD, NkyLlB, NkyLuB; 


extern int    NxLD, NyLD, NkyLD, NzLD, NvLD, NmLD, NsLD;
extern int    NxLB, NyLB, NzLB, NvLB, NmLB, NsLB;
extern int    NxGB, NyGB, NkyGD, NzGB, NvGB, NmGB, NsGB;

extern double dx, dy, dz, dv;
extern double Lx, Ly, Lz, Lv, Lm;

extern double *X, *Z, *V, *M, *Z;

extern int control_triggered_signal;

/****************************************************************************/



class FileIO;
class Timing;

class IfaceGKC {
   
    protected:
      virtual void printOn(std::ostream &output) const = 0;
      virtual ~IfaceGKC() { };
    public:
    friend std::ostream& operator<<(std::ostream& output, const IfaceGKC& ih) { ih.printOn(output); return output; };

    // Data Output Operation
};
   
// Forward declaration
class Plasma;
extern Plasma *plasma;

class Species;
extern Species *species;


typedef CComplex(*A6zz)[0][0][0][0][0];
typedef CComplex(*A5zz)[0][0][0][0];
typedef CComplex(*A4zz)[0][0][0];
typedef CComplex(*A3zz)[0][0];
typedef CComplex(*A2zz)[0];

typedef Real(*A2rr)[0];
typedef Real(*A3rr)[0][0];

#endif // __GLOBAL_H
