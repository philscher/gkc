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

#include<string>
#include <sstream>
#include <stddef.h>
#include "hdf5.h"
#include "hdf5_hl.h"

#include "FileIO.h"

#include "Plasma.h"

#include "Special/Vector3D.h"

typedef struct Complex_t {
      double r;   // real part
      double i;   // imaginary part
};
//}  Complex;



  // Constructor
  FileIO::FileIO(Parallel *_parallel, Setup *setup):     parallel(_parallel)
  {

    // Set Initial values
    inputFileName         = setup->get("DataOutput.InputFileName", "--- None ---");
    outputFileName        = setup->get("DataOutput.OutputFileName", "default.h5");
    info                  = setup->get("DataOutput.Info", "No information provided");
    resumeFile            = setup->get("DataOutput.Resume", 0);
    overwriteFile         = setup->get("DataOutput.Overwrite", 0) || (setup->flags & Setup::GKC_OVERWRITE);
   
     dataFileFlushTiming  = Timing(setup->get("DataOutput.Flush.Step", -1)       , setup->get("DataOutput.Flush.Time", 100.)); 
    

    ///////// Define Timeing Datatype
    timing_tid = H5Tcreate(H5T_COMPOUND, sizeof(Timing));
    H5Tinsert(timing_tid, "Timestep", HOFFSET(Timing, step), H5T_NATIVE_INT   );
    H5Tinsert(timing_tid, "Time"    , HOFFSET(Timing, time), H5T_NATIVE_DOUBLE);

    // don't changes r and i name otherwise it will break compatibiltiy with pyTables
    complex_tid = H5Tcreate(H5T_COMPOUND, sizeof (Complex_t));
    H5Tinsert(complex_tid, "r", HOFFSET(Complex_t,r), H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_tid, "i", HOFFSET(Complex_t,i), H5T_NATIVE_DOUBLE);
    
    vector3D_tid = H5Tcreate(H5T_COMPOUND, sizeof (Vector3D));
    H5Tinsert(vector3D_tid, "x", HOFFSET(Vector3D,x), H5T_NATIVE_DOUBLE);
    H5Tinsert(vector3D_tid, "y", HOFFSET(Vector3D,y), H5T_NATIVE_DOUBLE);
    H5Tinsert(vector3D_tid, "z", HOFFSET(Vector3D,z), H5T_NATIVE_DOUBLE);
 
    hsize_t species_dim[] = { setup->get("Grid.Ns", 1)}; 
    species_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, species_dim);
    
    // used for species name
          // BUG : Somehow HDF-8 stores up to 8 chars fro 64 possible. Rest are truncated ! Why ?
    s256_tid = H5Tcopy(H5T_C_S1); H5Tset_size(s256_tid, 64); H5Tset_strpad(s256_tid, H5T_STR_NULLTERM);


    offset0[0] = 0;
    offset0[1] = 0;
    offset0[2] = 0;
    offset0[3] = 0;
    offset0[4] = 0;
    offset0[5] = 0;
    offset0[6] = 0;
    offset0[7] = 0;
    

    // Create/Load HDF5 file
    if(resumeFile == false || (inputFileName != outputFileName)) create(setup);

   }

    int FileIO::create(Setup *setup) {
     
        hid_t file_plist = H5Pcreate(H5P_FILE_ACCESS);
#ifdef GKC_PARALLEL_MPI
   //       pass some information onto the underlying MPI_File_open call 
          MPI_Info file_info;
          check(MPI_Info_create(&file_info), DMESG("File info"));
          /* 
          H5Pset_sieve_buf_size(file_plist, 262144); 
          H5Pset_alignment(file_plist, 524288, 262144);
                
          MPI_Info_set(file_info, (char *) "access_style"        , (char *) "write_once");
          MPI_Info_set(file_info, (char *) "collective_buffering", (char *) "true");
          MPI_Info_set(file_info, (char *) "cb_block_size"       , (char *) "1048576");
          MPI_Info_set(file_info, (char *) "cb_buffer_size"      , (char *) "4194304");
           * */

          check( H5Pset_fapl_mpio(file_plist, parallel->Comm[DIR_ALL], file_info), DMESG("Set MPI Property"));
#endif
        file = check(H5Fcreate(outputFileName.c_str(), (overwriteFile ? H5F_ACC_TRUNC : H5F_ACC_EXCL),
                        H5P_DEFAULT, file_plist ), DMESG("H5FCreate : HDF5 File (File already exists ? use -f to overwrite) : " + outputFileName));
        check( H5Pclose(file_plist),   DMESG("H5Pclose"));

#ifdef GKC_PARALLEL_MPI
        MPI_Info_free(&file_info);
#endif
        
         //////////////////////////////////////////////////////////////// Info Group ////////////////////////////////////////////////////////

          hid_t infoGroup = check(H5Gcreate(file, "/Info",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));

         check(H5LTset_attribute_string(infoGroup, ".", "Output", outputFileName.c_str()), DMESG("H5LTset_attribute"));
         check(H5LTset_attribute_string(infoGroup, ".", "Input",  inputFileName.c_str()), DMESG("H5LTset_attribute"));
         
         
         check(H5LTset_attribute_string(infoGroup, ".", "Version", PACKAGE_VERSION), DMESG("H5LTset_attribute"));
         // Some Simulation specific stuff
         //check(H5LTset_attribute_string(infoGroup, ".", "Solver", ((setup->Solver & VL_LIN) ? "Linear" : "Non-Linear")), DMESG("H5LTset_attribute"));
         //heck(H5LTset_attribute_string(infoGroup, ".", "Type",   ((setup->VlasovType   & VLASOV_LOCAL ) ? "Local"  : "Global"    )), DMESG("H5LTset_attribute"));
         //heck(H5LTset_attribute_string(infoGroup, ".", "FFTSolverS",   ((setup->VlasovType   & VLASOV_LOCAL ) ? "Local"  : "Global"    )), DMESG("H5LTset_attribute"));
         //check(H5LTset_attribute_string(infoGroup, ".", "Initial Condition", setup->PerturbationMethod.c_str()), DMESG("H5LTset_attribute"));
         check(H5LTset_attribute_string(infoGroup, ".", "Info", info.c_str()), DMESG("H5LTset_attribute"));
         
         check(H5LTset_attribute_string(infoGroup, ".", "Config", setup->configFileString.c_str()), DMESG("H5LTset_attribute"));

         H5Gclose(infoGroup);
         
         
         /// Wrote setup constants, ugly here ////
         hid_t constantsGroup = check(H5Gcreate(file, "/Constants",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
         //
         if (!setup->parser_constants.empty()) { 
            
           std::vector<std::string> const_vec = Setup::split(setup->parser_constants, ",");

            for(int s = 0; s < const_vec.size(); s++) { 
                std::vector<std::string> key_value = Setup::split(const_vec[s],"=");
                double value = Setup::string_to_double(key_value[1]);
                int dim[] = { 1 };
   //           check(H5LTmake_dataset_double(constantsGroup, Setup::trimLower(key_value[0], false).c_str(), 1, dim, &value ), DMESG("Write Constants Attributes"));
                check(H5LTset_attribute_double(constantsGroup, ".", Setup::trimLower(key_value[0], false).c_str(), &value, 1), DMESG("H5LTset_attribute"));
                //check(H5LTset_attribute_double(constantsGroup, ".", Setup::trimLower(key_value[0], false).c_str(), &(Setup::string_to_double(key_value[1])), 1), DMESG("H5LTset_attribute"));
            };
         
         }
         
          H5Gclose(constantsGroup);

         
         // ********************* setup Table for CFL   *****************88
         cfl_table = new CFLTable();
         
         cfl_offset[0] =  HOFFSET( CFLTable, timeStep );
         cfl_offset[1] =  HOFFSET( CFLTable, time );
         cfl_offset[2] =  HOFFSET( CFLTable, Fx );
         cfl_offset[3] =  HOFFSET( CFLTable, Fy );
         cfl_offset[4] =  HOFFSET( CFLTable, Fz  );
         cfl_offset[5] =  HOFFSET( CFLTable, Fv );
         cfl_offset[6] =  HOFFSET( CFLTable, total );
          

         for(int i = 1; i < 7; i++)  cfl_sizes[i] = sizeof(double); cfl_sizes[0] = sizeof(int);
         hid_t   cfl_type[7]; for(int i = 1; i < 7; i++)  cfl_type [i] = H5T_NATIVE_DOUBLE; cfl_type[0] = H5T_NATIVE_INT;

         const char *cfl_names[7];
         cfl_names[0] = "timeStep";
         cfl_names[1] = "time";
         cfl_names[2] = "Fx"; cfl_names[3] = "Fy"; cfl_names[4] = "Fz"; cfl_names[5] = "Fv"; cfl_names[6] = "Total";

          check(H5TBmake_table("cflTable", file, "cfl", (hsize_t) 7, (hsize_t) 0, sizeof(CFLTable), (const char**) cfl_names,
                               cfl_offset, cfl_type, 32, NULL, 0, cfl_table ), DMESG("H5Tmake_table : cfl"));
         

         return GKC_SUCCESS;
    }


//    int FileIO::writeInitialConditions(Vlasov *vlasov, Fields *fields, Timing timing) {
    /* 
         hid_t initialGroup = check(H5Gcreate(file, "/Initial",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));

         hsize_t hs_Nx[1] = {Nx}, hs_Nv[1]={Nv};
         H5LTmake_dataset_double(initialGroup, "T" , 1, hs_Nx, T (RxLD).data());
         H5LTmake_dataset_double(initialGroup, "Tx", 1, hs_Nx, Tx(RxLD).data());
         H5LTmake_dataset_double(initialGroup, "N" , 1, hs_Nx, N (RxLD).data());

         H5Gclose(initialGroup);
     * */      
          
         // sometimes we want to skip IC to reduce file size

         

         
 //         return GKC_SUCCESS;
 //   }




      // This is a bit ugly ... :(
/* 
      int FileIO::load(Vlasov *vlasov, Fields *fields) {
        return GKC_FAILED;
        hid_t file_in;
        if(inputFileName == outputFileName)
           file_in = check(H5Fopen( inputFileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), DMESG("H5Fopen : inputFileName"));
    else

           file_in = check(H5Fopen( inputFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), DMESG("H5Fopen : inputFileName"));
          // E-Potential
          hid_t phi_plist_in = H5P_DEFAULT; 
#ifdef GKC_PARALLEL_MPI
          phi_plist_in = H5Pcreate(H5P_DATASET_XFER);
          H5Pset_dxpl_mpio(phi_plist_in, H5FD_MPIO_COLLECTIVE);
#endif     
          hid_t phi_dset_in = check(H5Dopen(file_in, "Phi", H5P_DEFAULT), DMESG(" H5Dopen : Phi"));
          H5LTget_attribute_int(phi_dset_in, ".", "phiCount", &phiCount);
          hid_t phi_dspace_in = H5Dget_space (phi_dset_in);
          // get last entry
          hsize_t phi_offset[] = {NzLlB-1, NyLlB-1, NxLlB-1, phiCount-1 }; 
        check(H5Sselect_hyperslab(phi_dspace_in, H5S_SELECT_SET, phi_offset, NULL, phi_chunkdim, NULL), DMESG("Selecting Hyperslab for phi"));
        check(H5Dread(phi_dset_in, H5T_NATIVE_DOUBLE, phi_mspace, phi_dspace_in, phi_plist, fields->phi.data()), DMESG("H5Dread : Phi"));
       
        // Phasespace function
        hid_t psf_plist_in = H5P_DEFAULT;
#ifdef GKC_PARALLEL_MPI
        psf_plist_in = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(psf_plist_in, H5FD_MPIO_COLLECTIVE);
#endif
        hid_t psf_dset_in = check(H5Dopen(file_in, "Phasespace", H5P_DEFAULT), DMESG(" H5Dopen : PSF"));
        check(H5LTget_attribute_int(psf_dset_in, ".", "psfCount", &psfCount), DMESG("H5LT_getattribute"));
        hid_t psf_dspace_in = H5Dget_space (psf_dset_in);
        hsize_t psf_offset[7] =  { 0, NmLlB-1,  NvLlB-1, NzLlB-1, NyLlB-1, NxLlB-1, psfCount-1 };
        check ( H5Sselect_hyperslab (psf_dspace_in, H5S_SELECT_SET, psf_offset, NULL, psf_chunkdim, NULL), DMESG("Selecting Hyperslab for PSF"));
        check(H5Dread(psf_dset_in, H5T_NATIVE_DOUBLE, psf_mspace, psf_dspace_in, psf_plist_in, vlasov->f.data()), DMESG("H5Dread : PSF"));
        
        // Get f0 ! .. definition the t=0 value is the f0 value !
        //psf_offset =  { NvLlB-1, NzLlB-1, NyLlB-1, NxLlB-1, 0 };
        psf_offset[6] = 0;
        check ( H5Sselect_hyperslab (psf_dspace_in, H5S_SELECT_SET, psf_offset, NULL, psf_chunkdim, NULL), DMESG("Selecting Hyperslab for PSF"));
        check(H5Dread(psf_dset_in, H5T_NATIVE_DOUBLE, psf_mspace, psf_dspace_in, psf_plist_in, vlasov->f0.data()), DMESG("H5Dread : PSF"));
        
       if (setup->inputFileName == setup->outputFileName) { file = file_in; }
        else {
             check( H5Fclose(file_in) , DMESG("Unable to close file ... exiting"));
             // reset position counters for new file
             phiCount=1; psfCount=1;
        }
            //file = check(H5Fopen( setup->outputFileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), DMESG("H5Fopen"));
           return GKC_SUCCESS;
    }

     

  int FileIO::writeCFLValues(Analysis *analysis, Fields *fields, Timing timing) {
      return GKC_SUCCESS;
        
    cfl_table->timeStep  = timing.step;
    cfl_table->time  = timing.time;
    cfl_table->Fx    = analysis->getMaxTimeStep(timing, DIR_X);
    cfl_table->Fy    = analysis->getMaxTimeStep(timing, DIR_Y);
    cfl_table->Fz    = analysis->getMaxTimeStep(timing, DIR_Z);
    cfl_table->Fv    = analysis->getMaxTimeStep(timing, DIR_V);
    cfl_table->total = analysis->getMaxTimeStep(timing);
    
    check(H5TBappend_records (file, "cfl", 1, sizeof(CFLTable), cfl_offset, cfl_sizes, cfl_table), DMESG("Append Table")); 

      return GKC_SUCCESS;
  }

  void zero(double *DArray, int N) {
      for(int i = 0; i < N; i++) DArray[i] = 0.e0;
  }

 * */








    // Destructor
    FileIO::~FileIO()  {
       // Free all HDF5 resources

     // close some extra stuff
     check( H5Tclose(complex_tid),   DMESG("H5Tclose"));
     check( H5Tclose(timing_tid),   DMESG("H5Tclose"));
     check( H5Tclose(species_tid),  DMESG("H5Tclose"));
     check( H5Tclose(s256_tid),     DMESG("H5Tclose"));
    
     
     // close file
     check( H5Fclose(file) , DMESG("Unable to close file ..."));

     delete cfl_table;

}
/* 
int FileIO::checkOutput(Vlasov *vlasov, Fields *fields, Visualization *visual, Analysis *analysis, Timing timing, double dt, bool force_write) {
        
        //if (timing.check(dataOutputPSF, dt) || force_write ) writePhaseSpace(vlasov->f, timing);
        if (timing.check(dataOutputVisual,dt) )              visual->doPlots(vlasov, fields, analysis, timing);
       if (timing.check(dataOutputStatistics, dt))          analysis->writeScalarValues(vlasov, fields,  analysis, timing);
       if (timing.check(dataOutputStatistics, dt))          writeCFLValues(analysis, fields, timing);


 * */ 

double FileIO::getOutputMaxTimeStep(Timing timing, double dt) {
       
    check(-1, DMESG("BLA"));        
  /* 
        dt = min(dt, dataOutputPSF        & timing);
        dt = min(dt, dataOutputPhi        & timing);
        dt = min(dt, dataOutputStatistics & timing);
        dt = min(dt, dataOutputXProperty  & timing);
   * */ 
//       if (timing % dataOutputStatistics) writeScalarValues(vlasov, fields,  analysis, timing);
//       if (timing % dataOutputXProperty ) writeXProperty(vlasov, fields,  analysis, timing);
//       if (timing % dataOutputStatistics) writeCFLValues(analysis, fields, timing);

        return dt; 

}
    
hid_t  FileIO::newGroup(std::string name, hid_t parentNode)
{
      if (parentNode == -2) parentNode = getFileID();
      hid_t newGroup = check(H5Gcreate(parentNode, name.c_str(),H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Spectrum : H5Gcreate"));
      return newGroup;

};


void FileIO::flush(Timing timing, double dt)
{
   if(timing.check(dataFileFlushTiming, dt)) H5Fflush(file, H5F_SCOPE_GLOBAL);

}


void FileIO::printOn(std::ostream &output) const {
         output << "            -------------------------------------------------------------------" << std::endl
                << "Data       |  Input     : " << inputFileName        << " Output    : " <<  outputFileName        << " Resume : " << ((resumeFile)?"yes":"no") << std::endl;
      };

FileAttr* FileIO::newTiming(hid_t group, hsize_t offset, bool write)
{
    
       hsize_t timing_chunkdim[] = {1}; hsize_t timing_maxdim[] = {H5S_UNLIMITED};
       hsize_t time_dim[1] = { 1 };
       FileAttr *Attr =  new FileAttr("Time", group, 1, time_dim, timing_maxdim, timing_chunkdim, &offset,  timing_chunkdim, &offset, parallel->myRank == 0, timing_tid);
    
       return Attr;

}
