/*
 * =====================================================================================
 *
 *       Filename: FileIO.h
 *
 *    Description: Data Input/Output using HDF-5 initialization routines
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FILEIO_H_
#define __FILEIO_H_

#include<fstream>
#include<iostream>
#include <hdf5.h>
#include <string>
#include <fstream>


#include "Parallel.h"
#include "Setup.h"
#include "Global.h"

#include "Timing.h"
#include "SHDF5/FileAttr.h"
#include "SHDF5/TableAttr.h"


class Visualization;

/**
*   @brief class for reading/writing data using HDF5
*
*  This class performs readind(resume previous simulation) and
*  writing Data. Additionally a Table is provided in which
*  Take care, we use Fortran Array syntax
*  This example writes data to the HDF5 file.
*  Data conversion is performed during write operation.  
*  
**/
class FileIO : public IfaceGKC 
{

  private:
   
   Timing dataFileFlushTiming; ///< Timing when to flush the HDF-5 file


   /**
   *  @brief Table to store the CFL values and various contributions
   *
   *  @note  Is it worth the effort to keep in the codebase ?
   *         Get rid of contributional part, just store the time-step
   **/ 
   typedef struct CFLTable
   {
     //  CFLTable(int _timeStep, double _time, double _Fx, double _Fy, double _Fz, double _total) : 
     //           timeStep(_timeStep), time(_time), Fx(_Fx), Fy(_Fy), Fz(_Fz), total(_total) {}; 
     int timeStep;
     double time;
     double   Fx, Fy, Fz, Fv, total;
   } _CFLTable;


   
   // Data files
   
   string inputFileName;     ///< Input file name (not working)
   string outputFileName;    ///< Name of the output file
   string info;              ///< additional information to append to the file
//   ofstream asciiOutputFile; 
   
   Parallel *parallel;

    
   size_t cfl_offset[7];
   size_t cfl_sizes[7];

   CFLTable *cfl_table;



public:
   /** Flush all data to disk to prevent corruption */
   void flush(Timing timing, double dt);
   
   hid_t timing_mspace, species_tid;
   
   hid_t s256_tid;
   
   hid_t file;  
   
   hsize_t offset0[10];//= { 0, 0, 0, 0, 0, 0, 0 };
   hid_t timing_tid, complex_tid, vector3D_tid;
   hid_t getFileID() const { return file; };

   // move to private
   bool resumeFile, overwriteFile;
   FileIO(Parallel *parallel, Setup *setup);
   virtual ~FileIO();

   //    int load(Vlasov *vlasov, Fields *fields);

   double getOutputMaxTimeStep(Timing timing, double dt);
   hid_t  newGroup(hid_t parentNode, std::string name);
        
   FileAttr *newTiming(hid_t group, hsize_t offset=0, bool write=1);

  protected:
   virtual void printOn(std::ostream &output) const;
   virtual void initDataOutput(FileIO *fileIO) {};
   virtual void writeData(Timing *timing) {};
   virtual void closeData() {};

  private:

   int create(Setup *setup);
};

#endif // __FILEIO_H_
