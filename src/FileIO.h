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

//! FileIO - class for reading/writing data using HDF5
/*!
 *  This class performs readind(resume previous simulation) and
 *  writing Data. Additionally a Table is provided in which
 *  Take care, we use Fortran Array syntax
   *  This example writes data to the HDF5 file.
   *  Data conversion is performed during write operation.  
 *  
 */
class FileIO : public IfaceHelios 
{
private:
Timing dataFileFlushTiming;

typedef struct CFLTable
{
//  CFLTable(int _timeStep, double _time, double _Fx, double _Fy, double _Fz, double _total) : 
//            timeStep(_timeStep), time(_time), Fx(_Fx), Fy(_Fy), Fz(_Fz), total(_total) {}; 
  int timeStep;
  double time;
  double   Fx, Fy, Fz, Fv, total;
} _CFLTable;


    // Data files
    string  inputFileName;
    string outputFileName;
    string info;
    // Timesteps to output the data 
//    Timing dataOutputPSF, dataOutputPhi, dataOutputVisual, dataOutputStatistics, dataOutputXProperty, dataOutputMoments;
    ofstream asciiOutputFile;

    Parallel *parallel;

     size_t cfl_offset[7];
     size_t cfl_sizes[7];

     CFLTable *cfl_table;





  protected:

     virtual void printOn(ostream &output) const {
         output << "            -------------------------------------------------------------------" << std::endl
                << "Data       |  Input     : " << inputFileName        << " Output    : " <<  outputFileName        << " Resume : " << ((resumeFile)?"yes":"no") << std::endl;
         /* 
                << "Output     | PSF        : " << dataOutputPSF       
                << "           | Phi        : " << dataOutputPhi      
                << "           | Statistics : " << dataOutputStatistics 
                << "           | XProperty  : " <<  dataOutputXProperty  
                << "           | Moments    : " << dataOutputMoments   
                << "           | Info       : " << info << std::endl
                << "           -------------------------------------------------------------------" << std::endl;
          * */
      };


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

//    int checkOutput(Vlasov *vlasov, Fields *fields, Visualization *visual, Analysis *analysis, Timing timing, double dt, bool force_write=false);
    double getOutputMaxTimeStep(Timing timing, double dt);
    hid_t  newGroup(hid_t parentNode, std::string name);
        
    FileAttr *newTiming(hid_t group, hsize_t offset=0, bool write=1) {
     hsize_t timing_chunkdim[] = {1}; hsize_t timing_maxdim[] = {H5S_UNLIMITED};
     hsize_t time_dim[1] = { 1 };
     FileAttr *Attr =  new FileAttr("Time", group, 1, time_dim, timing_maxdim, timing_chunkdim, &offset,  timing_chunkdim, &offset, parallel->myRank == 0, timing_tid);
     return Attr;
    }
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};
private:

    int create(Setup *setup);
};

#endif // __FILEIO_H_
