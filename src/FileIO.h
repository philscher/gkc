/*
 * =====================================================================================
 *
 *       Filename: FileIO.h
 *
 *    Description: Data Input/Output using HDF-5 (http://www.hdfgroup.org/HDF5/) 
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FILEIO_H_
#define __FILEIO_H_

#include "Global.h"

#include "Parallel/Parallel.h"
#include "Setup.h"
#include "Timing.h"

#include "SHDF5/FileAttr.h"
#include "SHDF5/TableAttr.h"


/**
*   @brief class for reading/writing data using HDF5
*
*  This class performs reading (resume previous simulation) and
*  writing Data. Additionally a Table is provided in which
*  Take care, we use Fortran Array syntax
*  This example writes data to the HDF5 file.
*  Data conversion is performed during write operation.  
*  
**/
class FileIO : public IfaceGKC 
{
  Parallel *parallel;

  Timing dataFileFlushTiming; ///< Timing when to flush the HDF-5 file
   
  std::string outputFileName,    ///< Name of the output file
              info;              ///< additional information to append to the file
  
  void create(Setup *setup, bool allowOverwrite);
  
 public:
  
  /**
  *   @brief constructor
  *
  *   Accepts following Setup parameters
  *
  *
  **/
  FileIO(Parallel *parallel, Setup *setup);
  
  virtual ~FileIO();
   
  std::string inputFileName;     ///< Input file name (move to protected!)

  /** 
  *
  * @brief Flush all data to disk to prevent corruption 
  *
  * Should be called in regular intervals in order to access written data
  * in case of program abort. This operation is of the order of ~seconds
  * thus should be used to often. 
  *
  **/
  void flush(Timing timing, double dt, bool force_flush=false);
   
  hid_t   complex_tid, ///< Complex Data type
           timing_tid, ///< Type id for Timing 
          species_tid, ///< species HDF-5 type id (where is it used ?)
        specfield_tid, ///< [ Species, 3 ] type for heat/particle flux
             str_tid,  ///< string data type
         vector3D_tid; ///< Vector (x,y,z) type
   
  hid_t file; ///< main data file id 
   
  hid_t getFileID() const { return file; };

  bool resumeFile; ///< true if simulation resumed from input file 

  /**
  *   @brief Created and returns a new group
  *
  *
  **/
  hid_t  newGroup(std::string name, hid_t parentNode=-2);
        
  /**
  *   @brief Returns a timing identifier
  *
  *
  **/
  FileAttr *newTiming(hid_t group, hsize_t offset=0, bool write=1);

 protected:

  virtual void printOn(std::ostream &output) const;
  virtual void initData(Setup *setup);
  virtual void writeData(const Timing &timing, const double dt) {};
  virtual void closeData() {};
};

#endif // __FILEIO_H_
