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


#define HELIOS_SUCCESS  1
#define HELIOS_FAILED   -1
#define HELIOS_STOP     2
#define HELIOS_WRITE     3
#define HELIOS_EXIT     4
#define HELIOS_FINISH     5


#if defined(__cplusplus)


#include<sstream>
#include<string>
#include<cstdlib>

#include<blitz/array.h>
#include<blitz/types.h>
#include<blitz/allocate.h>

/**
 *   Implicit OpenMP parallelization. Take care, it does not 
 *   handle any private or shared variables , or reduces or ...!!
 *
 *   HANDLE WITH CARE !
 * */
#define omp_for _Pragma("omp parallel for") for

using namespace blitz;

// DIR_SIZE has to be the number of elements, excluding DIR_SIZE
enum Direction : int {DIR_X=0, DIR_Y, DIR_Z, DIR_V, DIR_M, DIR_S, DIR_ALL, DIR_XYZ, DIR_VMS, DIR_MS, DIR_FFT, DIR_YZVMS, DIR_VM, DIR_YZ, DIR_XYZVM, DIR_XY, DIR_XMS, DIR_XM, DIR_SIZE};

      
    inline int check( int status, std::string file, int line, std::string error_text, bool segfault=false) {
        if(status == -1 ) {
            // check rank from MPI!
            std::stringstream ss;
            ss << std::endl; 
            ss << "\033[0;m";
            ss << "\033[1;m"  <<"!.... " << file << "(" << line << ") " << error_text; 
            ss << "\033[0m" << std::endl << std::flush;

            std::cout << ss.str();
              
            //TerminalIO::writeLogMessage(ss.str().c_str());
            // we produce SIGFPE (floating point exception), so we can backtrace stack
            if(segfault == true) {
                int *i = 0; *i = 0;
            }
            
	    abort();
            // exit through abort so we can get stack trace
#ifndef DEBUG
            exit(0);
#else
            exit(0);
#endif
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
extern Array1d X, Y, Z, V, M, kY;

// define some files for dataoutput
int writeMessage(const char *message);
int writeMessage(std::string message);




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

class IfaceHelios {

    protected:
      virtual void printOn(ostream &output) const = 0;
      //virtual ~IfaceHelios() { closeData(); };
      virtual ~IfaceHelios() { };
    public:
    friend ostream& operator<<(ostream& output, const IfaceHelios& ih) { ih.printOn(output); return output; };

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
extern GeneralArrayStorage<6> HeliosStorage;
extern GeneralArrayStorage<4> HeliosStorage4;


//#define HELIOS_GEOMETRY GeometrySlab
#define HELIOS_GEOMETRY Geometry2D

class HELIOS_GEOMETRY;



#endif // __GLOBAL_H
