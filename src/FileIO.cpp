/*
 * =====================================================================================
 *
 *       Filename: FileIO.cpp
 *
 *    Description: Data Input/Output using HDF-5 initialization routines
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include <string>
#include <sstream>
#include <stddef.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "FileIO.h"
#include "Plasma.h"
#include "Special/Vector3D.h"
#include "Tools/System.h"

typedef struct ComplexSplit_t {
  double r;   ///< real part
  double i;   ///< imaginary part
} ComplexSplit;


FileIO::FileIO(Parallel *_parallel, Setup *setup)  :  parallel(_parallel)
{

  // Set Initial values
  inputFileName        = setup->get("DataOutput.InputFileName" , "");
  outputFileName       = setup->get("DataOutput.OutputFileName", "default.h5");
  info                 = setup->get("DataOutput.Info"          , "No information provided");

  dataFileFlushTiming  = Timing(setup->get("DataOutput.Flush.Step", -1), setup->get("DataOutput.Flush.Time", 100.)); 
    
  resumeFile           = inputFileName != "";
  
  bool allowOverwrite  = setup->get("DataOutput.Overwrite", 0) || (setup->flags & Setup::GKC_OVERWRITE);
    

  ///////// Define Timing Datatype
  timing_tid = H5Tcreate(H5T_COMPOUND, sizeof(Timing));
  H5Tinsert(timing_tid, "Timestep", HOFFSET(Timing, step), H5T_NATIVE_INT   );
  H5Tinsert(timing_tid, "Time"    , HOFFSET(Timing, time), H5T_NATIVE_DOUBLE);

  // do not changes r and i name otherwise it will break compatibiltiy with pyTables
  complex_tid = H5Tcreate(H5T_COMPOUND, sizeof (ComplexSplit_t));
  H5Tinsert(complex_tid, "r", HOFFSET(ComplexSplit_t, r), H5T_NATIVE_DOUBLE);
  H5Tinsert(complex_tid, "i", HOFFSET(ComplexSplit_t, i), H5T_NATIVE_DOUBLE);
    
  vector3D_tid = H5Tcreate(H5T_COMPOUND, sizeof (Vector3D));
  H5Tinsert(vector3D_tid, "x", HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
  H5Tinsert(vector3D_tid, "y", HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
  H5Tinsert(vector3D_tid, "z", HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);

  hsize_t species_dim[1]   = { setup->get("Grid.Ns", 1) }; 
  species_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, species_dim);
  
  hsize_t specfield_dim[2] = { setup->get("Grid.Ns", 1), setup->get("Plasma.Beta", 0.) == 0. ? 1 : 2 }; 
  specfield_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, specfield_dim);
    
  // used for species name
  // BUG : Somehow HDF-8 stores only up to 8 chars of 64 possible. The rest is truncated ! Why ?
  // not supported by HDF-5 parallel yet ....
  // str_tid = H5Tcopy(H5T_C_S1); H5Tset_size(str_tid, H5T_VARIABLE); H5Tset_strpad(str_tid, H5T_STR_NULLTERM);
  str_tid = H5Tcopy(H5T_C_S1); H5Tset_size(str_tid, 64); H5Tset_strpad(str_tid, H5T_STR_NULLTERM);

  // Create/Load HDF5 file
  if(resumeFile == false || (inputFileName != outputFileName)) create(setup, allowOverwrite);

}


void FileIO::create(Setup *setup, bool allowOverwrite) 
{

  hid_t file_apl = H5Pcreate(H5P_FILE_ACCESS);

  // from http://mail.hdfgroup.org/pipermail/hdf-forum_hdfgroup.org/2010-June/012076.html
  // to prevent corruption of HDF-5 file (requires regular calls to this->flush() ] 
  {
//    H5AC_cache_config_t mdc_config;

//    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
//    H5Pget_mdc_config(file_apl, &mdc_config);
  
//    mdc_config.evictions_enabled = false;
//    mdc_config.incr_mode         = H5C_incr__off;
//    mdc_config.decr_mode         = H5C_decr__off;

//    H5Pset_mdc_config(file_apl, &mdc_config);
  }
  
  // pass some information onto the underlying MPI_File_open call 
  // to optimize file access
  {
#ifdef GKC_PARALLEL_MPI
    MPI_Info file_info;
    check(MPI_Info_create(&file_info), DMESG("File info"));
     
    // H5Pset_sieve_buf_size(file_plist, 262144); 
    // H5Pset_alignment(file_plist, 524288, 262144);
     
    // MPI_Info_set(file_info, (char *) "access_style"        , (char *) "write_once");
    // MPI_Info_set(file_info, (char *) "collective_buffering", (char *) "true");
    // MPI_Info_set(file_info, (char *) "cb_block_size"       , (char *) "1048576");
    // MPI_Info_set(file_info, (char *) "cb_buffer_size"      , (char *) "4194304");

    check( H5Pset_fapl_mpio(file_apl, parallel->Comm[DIR_ALL], file_info), DMESG("Critical HDF-5 Error"));

    // shouldn't be freed before H5Pclose(file_apl)
    MPI_Info_free(&file_info);
#endif
  }

  // Close file even if some objects are still open (should generate warning as this is not correct
  H5Pset_fclose_degree(file_apl, H5F_CLOSE_STRONG);
 
  // Create new output file
  file = check(H5Fcreate(outputFileName.c_str(), (allowOverwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL),
          H5P_DEFAULT, file_apl ), DMESG("H5FCreate : HDF5 File (File already exists ? use -f to overwrite)"));
     
  check( H5Pclose(file_apl),   DMESG("H5Pclose"));

        
  //////////////////////////////////////////////////////////////// Info Group ////////////////////////////////////////////////////////
   
  hid_t infoGroup = newGroup("/Info");
   
  check(H5LTset_attribute_string(infoGroup, ".", "Output" , outputFileName.c_str()), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_string(infoGroup, ".", "Input"  , inputFileName.c_str()) , DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_string(infoGroup, ".", "Version", PACKAGE_VERSION)       , DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_string(infoGroup, ".", "Info"   , info.c_str())          , DMESG("H5LTset_attribute"));

  check(H5LTset_attribute_string(infoGroup, ".", "Config", setup->configFileString.c_str()), DMESG("H5LTset_attribute"));
  
  check(H5LTset_attribute_string(file, ".", "StartTime", System::getTimeString().c_str()), DMESG("H5LTset_attribute"));

  // get & save HDF-5 version
  { 
    std::stringstream version;
    unsigned int majnum, minnum, relnum;
    H5get_libversion(&majnum, &minnum, &relnum);
    version << majnum << "." << minnum << "." << relnum << std::endl;
  
    check(H5LTset_attribute_string(infoGroup, ".", "HDF5version", version.str().c_str()), DMESG("H5LTset_attribute"));
  }

  H5Gclose(infoGroup);
         
  ////////////// Wrote setup constants, ugly here /////////////////
  hid_t constantsGroup = newGroup("/Constants");

  // Should be in Setup (Store Setup Constants
  if (!setup->parser_constants.empty()) { 
   
    std::vector<std::string> const_vec = Setup::split(setup->parser_constants, ",");

    // for(auto key_value : Setup::split(setup->parser_constants, ","))
    for(int s = 0; s < const_vec.size(); s++) { 
     
       std::vector<std::string> key_value = Setup::split(const_vec[s],"=");

       double value = std::stod(key_value[1]);
       
       check(H5LTset_attribute_double(constantsGroup, ".", key_value[0].c_str(), &value, 1), DMESG("H5LTset_attribute"));
    }
         
  }
   
  H5Gclose(constantsGroup);

}


void FileIO::initData(Setup *setup)
{

}


// Destructor
FileIO::~FileIO()  
{
   
  check(H5LTset_attribute_string(file, ".", "StopTime", System::getTimeString().c_str()), DMESG("H5LTset_attribute"));
  // Free all HDF5 resources

  // close some extra stuff
  check( H5Tclose(complex_tid), DMESG("H5Tclose"));
  check( H5Tclose(timing_tid ), DMESG("H5Tclose"));
  check( H5Tclose(species_tid), DMESG("H5Tclose"));
  check( H5Tclose(str_tid    ), DMESG("H5Tclose"));

  // close file
  check( H5Fclose(file)    , DMESG("Unable to close file ..."));

}

hid_t  FileIO::newGroup(std::string name, hid_t parentNode)
{
   if (parentNode == -2) parentNode = getFileID();
   hid_t newGroup = check(H5Gcreate(parentNode, name.c_str(),H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), 
                          DMESG("Error creating group file for Spectrum : H5Gcreate"));
   return newGroup;
};


void FileIO::flush(Timing timing, double dt)
{
   if(timing.check(dataFileFlushTiming, dt)) H5Fflush(file, H5F_SCOPE_GLOBAL);
}


void FileIO::printOn(std::ostream &output) const 
{
   output << "            -------------------------------------------------------------------" << std::endl
          << "Data       |  Input : " << (inputFileName == "" ? "---None---" : inputFileName)  << std::endl 
          << "           | Output : " <<  outputFileName  << " Resume : " << (resumeFile ? "yes" : "no") << std::endl;
};

FileAttr* FileIO::newTiming(hid_t group, hsize_t offset, bool write)
{
    
   hsize_t timing_cdim  [1] = {1             };
   hsize_t timing_maxdim[1] = {H5S_UNLIMITED };
   hsize_t timing_dim   [1] = {1             };
   
   FileAttr *Attr = new FileAttr("Time", group, file, 1, timing_dim, timing_maxdim, timing_cdim, &offset,  
                                  timing_cdim, &offset, parallel->myRank == 0, timing_tid);
    
   return Attr;
}

