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
#define GKC_FAILED   -1
#define GKC_STOP     2
#define GKC_WRITE    3
#define GKC_EXIT     4
#define GKC_FINISH   5

#include <cilk/cilk.h>


// Needed in case we include into Fortran
#if defined(__cplusplus)


//#include<sstream>
#include<string>
//#include<cstdlib>


// Define Scalar types (only for doubles now)
//#define GKC_C_COMPLEX

#ifdef GKC_C_COMPLEX
#include <complex.h>
typedef _Complex double      Complex;  
typedef double               Real   ;
inline Real abs(Complex A) { return 1.; };
inline Complex exp(Complex A) { return 1.; };

//#define newComplex(A,B) A + B * I
//#define _Complex_I ((double _Complex){0, 1})
#define newComplex(A,B) (A + B)
//#define Complex(A,B) (A + B)
#else
#include <complex>
typedef std::complex<double> Complex;  
typedef double               Real   ;
#define newComplex(A,B) Complex(A,B)
#endif

#include<blitz/array.h>
#include<blitz/allocate.h>
#include<blitz/types.h>

typedef _Complex double CComplex;  

#define _Imaginary ((CComplex) (0.+1.j)); 
//typedef (__extension__ 1.0i) Imaginary; 


typedef blitz::Array<Complex, 1>  Array1C;
typedef blitz::Array<Complex, 2>  Array2C;
typedef blitz::Array<Complex, 4>  Array4C;
typedef blitz::Array<Complex, 3>  Array3C;
typedef blitz::Array<Complex, 5>  Array5C;
typedef blitz::Array<Complex, 6>  Array6C;

typedef blitz::Array<double, 1>  Array1R;
typedef blitz::Array<double, 2>  Array2R;
typedef blitz::Array<double, 3>  Array3R;
typedef blitz::Array<double, 4>  Array4R;
typedef blitz::Array<double, 5>  Array5R;
typedef blitz::Array<double, 6>  Array6R;

using namespace blitz;

/**
 *   Implicit OpenMP parallelization. Take care, it does not 
 *   handle any private or shared variables , or reduces or ...!!
 *
 *   HANDLE WITH CARE !
 **/
#define omp_for  _Pragma("omp parallel for") for
#define simd_for _Pragma("simd") for



// DIR_SIZE has to be the number of elements, excluding DIR_SIZE
//enum Direction : int {DIR_X=0, DIR_Y, DIR_Z, DIR_V, DIR_M, DIR_S, DIR_ALL, DIR_XYZ, DIR_VMS, DIR_MS, DIR_FFT, DIR_YZVMS, DIR_VM, DIR_YZ, DIR_XYZVM, DIR_XY, DIR_XMS, DIR_XM, DIR_SIZE};
///enum Direction : int {DIR_X=0, DIR_Y, DIR_Z, DIR_V, DIR_M, DIR_S, DIR_ALL, DIR_XYZ, DIR_VMS, DIR_MS, DIR_VM, DIR_XYZVM, DIR_XY, DIR_XMS, DIR_XM, DIR_SIZE};
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

extern int NkyLlD, NkyLuD;

extern Range RzLD, RyLD, RxLD, RvLD, RmLD, RsLD, RkyLD; 
extern Range RxLB, RyLB, RzLB, RvLB, RmLB, RsLB;
extern Range RxLB4, RyLB4; 
extern Range RB, RB4, RFields ; 
extern Array1R X, Y, Z, V, M;

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
extern double dx, dy, dz, dv, dm;


/****************************************************************************/


// BUG find a way to remove them


extern int process_rank;


// extern stuff for Fortran Interface
//
extern bool do_gyro;


template<typename T> std::string Num2String(T number) {
    std::stringstream ss;
    ss << number;
    return ss.str();
}

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
extern GeneralArrayStorage<6> GKCStorage;
extern GeneralArrayStorage<4> GKCStorage4;
extern GeneralArrayStorage<3> GKCStorage3;


typedef CComplex(*A6zz)[][][][][];
typedef CComplex(*A5zz)[][][][];
typedef CComplex(*A4zz)[][][];
typedef CComplex(*A3zz)[][];
typedef CComplex(*A2zz)[][];
typedef CComplex(*A1zz)[];

typedef Complex(*A6z)[][][][][];
typedef Complex(*A5z)[][][][];
typedef Complex(*A4z)[][][];
typedef Complex(*A3z)[][];
typedef Complex(*A2z)[][];






#endif // __GLOBAL_H
