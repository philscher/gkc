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


#define GKC_SUCCESS  1

#include <cilk/cilk.h>


#include<string>
#include <complex>
#include <iostream>
typedef std::complex<double> Complex;  
typedef double               Real   ;


#include<blitz/array.h>

#include "external/allocate.h"

typedef _Complex double CComplex;  

#define _Imaginary ((CComplex) (0.+1.j)); 

extern "C" double cabs  (CComplex z);
extern "C" double creal (CComplex z);
extern "C" double cimag (CComplex z);
extern "C" double carg  (CComplex z);
extern "C" CComplex conj(CComplex z);
extern "C" CComplex cexp (CComplex z);


using namespace std;

template<class T> inline T pow2(T x) { return x*x; };
template<class T> inline T pow3(T x) { return x*x*x; };
template<class T> inline T pow4(T x) { const T x2 =  (x*x); return x2*x2; };
template<class T> inline T pow5(T x) { const T x2 = x * x; return x2*x2*x; };
template<class T> inline T pow6(T x) { const T x2 = x * x; return x2*x2*x2; };
template<class T> inline T pow8(T x) { const T x2 = x * x; x2 *= x2; return x2*x2; };
//__declspec(vector) inline CComplex square(CComplex x) { return x*x; };
inline CComplex square(CComplex x) { return x*x; };



/* 
int ipow(int base, int exp)
{
   int result = 1;
   while (exp)
   {
      
     if (exp & 1)
     
       result *= base;
       
     exp >>= 1;
     
     
     base *= base;
     
   }

              return result;
}
 * */


/**
 *   Implicit OpenMP parallelization. Take care, it does not 
 *   handle any private or shared variables , or reduces or ...!!
 *
 *   HANDLE WITH CARE !
 **/
#define omp_for     _Pragma("omp parallel for") for
// warning This is an OpenMP 3.0 function, is it possible to check for it ? (_OPENMP)
#define omp_C2_for  _Pragma("omp parallel for collapse(2)") for
#define omp_C3_for  _Pragma("omp parallel for collapse(3)") for
// for vectorization support, see 
#define simd_for _Pragma("simd") for


// DIR_SIZE has to be the number of elements, excluding DIR_SIZE
enum Dir : int {DIR_X=0, DIR_Y, DIR_Z, DIR_V, DIR_M, DIR_S, DIR_ALL, DIR_XYZ, DIR_VMS, DIR_MS, DIR_VM, DIR_XY, DIR_XYZVM, DIR_SIZE};
     
inline int check( int status, std::string file, int line, std::string error_text, bool doAbort=false) {
        if(status == -1 ) {
            
       // check rank from MPI!
            std::stringstream ss;
            ss << std::endl; 
            ss << "\033[0;m";
            ss << "\033[1;m"  <<"!.... " << file << "(" << line << ") " << error_text; 
            ss << "\033[0m" << std::endl << std::flush;

            std::cout << ss.str();
              
            // exit through abort so we can get stack trace
       if(doAbort==true) abort();
            abort();
            exit(0);
        }
        return status;
    }

    inline bool check( bool status, std::string file, int line, std::string error_text) {
        if(status == false) check( -1, file, line, error_text);
        return status;
}



/// Deprecated
extern blitz::Range RzLD, RyLD, RxLD, RvLD, RmLD, RsLD, RkyLD; 
extern blitz::Range RxLB, RyLB, RzLB, RvLB, RmLB, RsLB;
extern blitz::Range RxLB4, RyLB4; 
extern blitz::Range RB, RB4,  RFields ; 


// use assert instead !
#define DMESG(mesg)  std::string(__FILE__), __LINE__, std::string(mesg)
#define _D(mesg, value) std::cout << mesg << " " << value << std::endl << std::flush


#include "config.h"


/*************************  Global Variables   ********************************/

// DO NOT CHANGE THESE VARIABLES BY ANY MEANS OUTSIDE GRID
// OR CONFIG !!

// Gx[LlD], Gx[LuD], Gx[LuB], Gx[LuD] // should improve speed due to caching

extern int NkyLlD, NkyLuD;

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


extern double Lx, Ly, Lz, Lv, Lm;
extern int    Nx, Nky, Nz, Nv, Nm, Ns;
extern int    NxLD, NyLD, NkyLD, NzLD, NvLD, NmLD, NsLD;
extern int    NxLB, NyLB, NkyLB, NzLB, NvLB, NmLB, NsLB;
extern int    NxGB, NyGB, NkyGB, NzGB, NvGB, NmGB, NsGB;
extern double dx, dy, dz, dv;

extern double *X, *Y, *Z, *V, *M, *Z;


/****************************************************************************/





#include "hdf5.h"
#include "hdf5_hl.h"


class FileIO;
class Timing;

class IfaceGKC {
   
    protected:
      virtual void printOn(std::ostream &output) const = 0;
      //virtual ~IfaceGKC() { closeData(); };
      virtual ~IfaceGKC() { };
    public:
    friend std::ostream& operator<<(std::ostream& output, const IfaceGKC& ih) { ih.printOn(output); return output; };

    // Data Output Operation
};
    
class IfaceDataIO {
    virtual ~IfaceDataIO() { };
    virtual void initDataOutput(FileIO *fileIO) = 0;
    virtual void writeData(Timing timing, double dt) = 0;
    virtual void closeData() = 0;
};

// Forward declaration
class Plasma;
extern Plasma *plasma;


typedef CComplex(*A6zz)[][][][][];
typedef CComplex(*A5zz)[][][][];
typedef CComplex(*A4zz)[][][];
typedef CComplex(*A3zz)[][];
typedef CComplex(*A2zz)[];

typedef Real(*A2rr)[];


#endif // __GLOBAL_H
