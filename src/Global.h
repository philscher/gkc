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
#define GKC_WRITE    3

#include <cilk/cilk.h>


// Needed in case we include into Fortran
#if defined(__cplusplus)

#include<string>
#include <complex>
typedef std::complex<double> Complex;  
typedef double               Real   ;
#define newComplex(A,B) Complex(A,B)


#include<blitz/array.h>
#include<blitz/allocate.h>
#include<blitz/types.h>

typedef _Complex double CComplex;  

#define _Imaginary ((CComplex) (0.+1.j)); 
//typedef (__extension__ 1.0i) Imaginary; 

extern "C" double cabs  (CComplex z);
extern "C" double creal (CComplex z);
extern "C" double cimag (CComplex z);
extern "C" double carg  (CComplex z);
extern "C" CComplex conj(CComplex z);
extern "C" CComplex cexp (CComplex z);


typedef blitz::Array<Complex, 1>  Array1C;
typedef blitz::Array<Complex, 2>  Array2C;
typedef blitz::Array<Complex, 4>  Array4C;
typedef blitz::Array<Complex, 3>  Array3C;
typedef blitz::Array<Complex, 5>  Array5C;
typedef blitz::Array<Complex, 6>  Array6C;

typedef blitz::Array<double, 1>  Array1R;
typedef blitz::Array<double, 2>  Array2R;

using namespace blitz;

//template<class T> inline T pow2(T x) { return x*x; };

/**
 *   Implicit OpenMP parallelization. Take care, it does not 
 *   handle any private or shared variables , or reduces or ...!!
 *
 *   HANDLE WITH CARE !
 **/
#define omp_for     _Pragma("omp parallel for") for
// warning This is an OpenMP 3.0 function, is it possible to check for it ? (_OPENMP)
#define omp_for_C2  _Pragma("omp parallel for collapse(2)") for
#define omp_for_C3  _Pragma("omp parallel for collapse(3)") for
// for vectorization support, see 
#define simd_for _Pragma("simd") for


// DIR_SIZE has to be the number of elements, excluding DIR_SIZE
enum Direction : int {DIR_X=0, DIR_Y, DIR_Z, DIR_V, DIR_M, DIR_S, DIR_ALL, DIR_XYZ, DIR_VMS, DIR_MS, DIR_VM, DIR_XY, DIR_XYZVM, DIR_SIZE};
     

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
extern blitz::Range RB, RB4, RFields ; 


// define some files for dataoutput
int writeMessage(const char *message);
int writeMessage(std::string message);



// use assert instead !
#define DMESG(mesg)  std::string(__FILE__), __LINE__, std::string(mesg)
#define _D(mesg, value) std::cout << mesg << " " << value << std::endl << std::flush


#include "config.h"


#endif // __cplusplus

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
      virtual void printOn(ostream &output) const = 0;
      //virtual ~IfaceGKC() { closeData(); };
      virtual ~IfaceGKC() { };
    public:
    friend ostream& operator<<(ostream& output, const IfaceGKC& ih) { ih.printOn(output); return output; };

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


// needed for signal handling
extern blitz::GeneralArrayStorage<6> GKCStorage;
extern blitz::GeneralArrayStorage<4> GKCStorage4;
extern blitz::GeneralArrayStorage<3> GKCStorage3;


typedef CComplex(*A6zz)[][][][][];
typedef CComplex(*A5zz)[][][][];
typedef CComplex(*A4zz)[][][];
typedef CComplex(*A3zz)[][];
typedef CComplex(*A2zz)[];
//typedef CComplex(*A1zz)[];

typedef Real(*A2rr)[][];

typedef Complex(*A6z)[][][][][];
typedef Complex(*A5z)[][][][];
typedef Complex(*A4z)[][][];
typedef Complex(*A3z)[][];
typedef Complex(*A2z)[][];







#endif // __GLOBAL_H
