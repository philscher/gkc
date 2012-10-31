/*
 * =====================================================================================
 *
 *       Filename:  DataOutput.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/10/2011 03:18:44 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "DataOutput.h"

#include <iostream>

struct Complex_t {
    double r;
    double i;
} Complex_H5;

extern "C" double cabs(const cmplxd a);


struct Timing_t {
    int Step;
    double Time;
} Timing_H5;
   
DataOutput::DataOutput(Input *input, int Nx, int Nky, double Lx, double Ly, double *X, double *Psi0, double Viscosity, double Resistivity) 
{


     OutputTimeField = input->get("OutputTimeField", 25.);
     OutputTimePower = input->get("OutputTimePower", 1.);


      std::string outputFileName(input->get("OutputFilename", "Output.h5"));
      fileOutput = check(H5Fcreate(outputFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT ), DMESG("H5FCreate : HDF5 File (File already exists ? use -f to overwrite) : " + outputFileName));
        
     // don't changes r and i name otherwise it will break compatibiltiy with pyTables
     hid_t complex_tid = H5Tcreate(H5T_COMPOUND, sizeof (Complex_H5));
     H5Tinsert(complex_tid, "r", HOFFSET(Complex_t,r), H5T_NATIVE_DOUBLE);
     H5Tinsert(complex_tid, "i", HOFFSET(Complex_t,i), H5T_NATIVE_DOUBLE);

      
       // Write Initial Conditions

       // Define ranges for normal outputs
       hsize_t x_dim[]       = { Nx, 1 };
       hsize_t x_maxdim[]    = { Nx, H5S_UNLIMITED };
       hsize_t x_chunkBdim[] = { Nx, 1 };
       hsize_t x_chunkdim[]  = { Nx, 1 };
       hsize_t x_moffset[]   = { 0 , 0 };
       
       // Define ranges for normal outputs
       hsize_t ky_dim[]       = { Nky, 1 };
       hsize_t ky_maxdim[]    = { Nky, H5S_UNLIMITED };
       hsize_t ky_chunkdim[]  = { Nky, 1 };
       hsize_t ky_chunkBdim[] = { Nky, 1 };
       hsize_t ky_moffset[]   = { 0  , 0 };

       // Define ranges for xy- 2d output (we also have complex conjugates)
       hsize_t xky_dim[]       = {   Nky, Nx,  1 };
       hsize_t xky_maxdim[]    = {   Nky, Nx, H5S_UNLIMITED };
       hsize_t xky_chunkBdim[] = {   Nky, Nx, 1 };
       hsize_t xky_chunkdim[]  = {   Nky, Nx, 1 };
       hsize_t xky_moffset[]   = {   0  ,  0, 0 };
       
      
       // For Time
       hid_t timing_tid = H5Tcreate(H5T_COMPOUND, sizeof(Timing_H5));
       H5Tinsert(timing_tid, "Step",     HOFFSET(Timing_t, Step), H5T_NATIVE_INT   );
       H5Tinsert(timing_tid, "Time"    , HOFFSET(Timing_t, Time), H5T_NATIVE_DOUBLE);
 
       hsize_t timing_chunkdim[] = {1}; hsize_t timing_maxdim[] = {H5S_UNLIMITED};
       hsize_t time_dim[1] = { 1 }; 
       hsize_t offset0 [] = { 0, 0, 0};

       hid_t initGroup = check(H5Gcreate(fileOutput, "/Init",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
       
       check(H5LTset_attribute_double(initGroup, ".", "X"  , X   , Nx), DMESG("Write Constants Attributes"));
       check(H5LTset_attribute_double(initGroup, ".", "Lx" , &Lx , 1 ), DMESG("Write Constants Attributes"));
       check(H5LTset_attribute_double(initGroup, ".", "Ly" , &Ly , 1 ), DMESG("Write Constants Attributes"));
       check(H5LTset_attribute_int   (initGroup, ".", "Nx" , &Nx , 1 ), DMESG("Write Constants Attributes"));
       check(H5LTset_attribute_int   (initGroup, ".", "Nky", &Nky, 1 ), DMESG("Write Constants Attributes"));

       check(H5LTset_attribute_double(initGroup, ".", "Viscosity", &Viscosity, 1 ), DMESG("Write Constants Attributes"));
       check(H5LTset_attribute_double(initGroup, ".", "Resistivity", &Resistivity, 1 ), DMESG("Write Constants Attributes"));
       H5Gclose(initGroup);
       
       
       hid_t fieldGroup = check(H5Gcreate(fileOutput, "/Field",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
       
       FA_Psi_Time  = new FileAttr("Time", fieldGroup, 1, time_dim, timing_maxdim, timing_chunkdim, offset0   , timing_chunkdim, offset0, 1, timing_tid); 
       FA_Psi       = new FileAttr("Psi" , fieldGroup, 3,  xky_dim,    xky_maxdim,    xky_chunkdim, xky_moffset,  xky_chunkBdim, offset0, 1, complex_tid);
       FA_Phi       = new FileAttr("Phi" , fieldGroup, 3,  xky_dim,    xky_maxdim,    xky_chunkdim, xky_moffset,  xky_chunkBdim, offset0, 1, complex_tid);
       
       FileAttr 
          *FA_Psi0 = new FileAttr("Psi0", fieldGroup, 2,    x_dim,      x_maxdim,      x_chunkdim,   x_moffset,    x_chunkBdim, offset0, 1); 
       FA_Psi0->write(Psi0);
       delete FA_Psi0;
       
       H5Gclose(fieldGroup);




       //  For Power-Spectra
       hid_t powerGroup = check(H5Gcreate(fileOutput, "/Power",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
       FA_PowerKy_Time  = new FileAttr("Time", powerGroup,  1, time_dim, timing_maxdim, timing_chunkdim, offset0, timing_chunkdim, offset0, 1, timing_tid); 
       FA_PowerKy_Psi   = new FileAttr("PsiY"   , powerGroup,  2,   ky_dim,     ky_maxdim,     ky_chunkdim, offset0,     ky_chunkdim, offset0, 1); 
       FA_PowerKy_Phi   = new FileAttr("PhiY"   , powerGroup,  2,   ky_dim,     ky_maxdim,     ky_chunkdim, offset0,     ky_chunkdim, offset0, 1); 
       
       H5Gclose(powerGroup);
     

       // Write

    };



DataOutput::~DataOutput() {

    delete FA_Psi_Time; delete FA_Psi, delete FA_Phi;
    delete FA_PowerKy_Time; delete FA_PowerKy_Phi; delete FA_PowerKy_Psi;
    H5Fclose(fileOutput);
   
}


void DataOutput::output(int Step, double Time, double dt, int Nky, int Nx, cmplxd Phi[Nky][Nx], cmplxd Psi[Nky][Nx]) {

       // write Fields         
      if(((fmod(Time+dt, OutputTimeField) < dt))) {
                    std::cout << " Writing Fields" << std::endl;
                    FA_Psi->write(&Psi[0][0]);
                    FA_Phi->write(&Phi[0][0]);
                    Timing_t T; T.Time = Time; T.Step = Step;
                    FA_Psi_Time->write(&T);
      }

      // write  Power
       // write Fields         
      if(((fmod(Time+dt, OutputTimePower) < dt))) {

         double PowerKy_Phi[Nky], PowerKy_Psi[Nky];
                    
         #pragma omp parallel for 
         for(int y_k = 0; y_k < Nky ; y_k++) { 
         
           PowerKy_Psi[y_k] = __sec_reduce_add(cabs(Psi[y_k][:]));
           PowerKy_Phi[y_k] = __sec_reduce_add(cabs(Phi[y_k][:]));

         }
                    
         FA_PowerKy_Psi->write(PowerKy_Psi); 
         FA_PowerKy_Phi->write(PowerKy_Phi); 
         Timing_t T; T.Time = Time; T.Step = Step;
         FA_PowerKy_Time->write(&T);

            
         std::cout << " Writing Power : psi : " << __sec_reduce_add(PowerKy_Psi[:])
                                   << " phi : " << PowerKy_Phi[1]
                                   //<< " phi : " << __sec_reduce_add(PowerKy_Phi[:])
                                   <<  std::endl;

      }


}

