/*
 * =====================================================================================
 *
 *       Filename: Grid.h
 *
 *    Description: Grid definitions incuding boundary
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GRID_H_
#define __GRID_H_


#include "config.h"
#include "Setup.h"
#include "Parallel.h"
#include "Integration.h"
#include "Global.h"

#include "FileIO.h"
//! Grid - Class setting basic properties and size of the computational domain
/*! The compuational domain consist basically of..
*   \image html grid.png
    For the Parallized these boudaries are transfered between the
*   CPUs using MPI. We use 2-D decomposition in X and Y. Because the
*
*
*/

class Grid : public IfaceHelios {

public:
Array4d  SendYu, SendXu, SendYl, SendXl; 
Array4d  RecvYu, RecvXu, RecvYl, RecvXl;

Integrate *Ipol_M;

// GhostCell
int NxGC, NyGC, NzGC, NvGC;

Range RxGB, RyGB, RzGB, RvGB, RmGB, RsGB; 
Range RxGD, RyGD, RzGD, RvGD, RmGD, RsGD; 


// Local Grid size (for Domain & Boundary)

      int NxGD, NyGD, NkyGD, NzGD, NvGD, NmGD, NsGD;


  Grid(Setup *setup, Parallel *parallel, FileIO *fileIO);
 ~Grid();
protected:
     virtual void printOn(ostream &output) const {
         output   << " Domain    |  Lx : " << Lx << "  Ly : " << Ly << "  Lz : " << Lz << "  Lv : " << Lv << ((do_gyro) ? std::string("  Lm : ") + Num2String(Lm) : "") << std::endl
                  << " Grid      |  Nx : " << Nx << "  Nky : " << Nky << "  Nz : " << Nz << "  Nv : " << Nv << ((do_gyro) ? std::string("  Nm : ") + Num2String(Nm) : "") << std::endl;
     
            
 #ifdef _DEBUG
            output << 
              "Domain     | " <<
            
                        "X(" << Nx << " " << NxGlD << "-" << NxGuD << ") "  
                        "Y(" << Ny << " " << NyGlD << "-" << NyGuD << ") "  
                        "Z(" << Nz << " " << NzGlD << "-" << NzGuD << ") "  
                        "V(" << Nv << " " << NvGlD << "-" << NvGuD << ") "  
                        "M(" << Nm << " " << NmGlD << "-" << NmGuD << ") "  
                        "S(" << Ns << " " << NsGlD << "-" << NsGuD << ") "  
                                                                            << std::endl
                        << "          | " << std::endl ;
#endif
           
    }
public: 
    void initDataOutput(FileIO *fileIO) {
          
      hid_t gridGroup = check(H5Gcreate(fileIO->getFileID(), "/Grid",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
          

         
          // Length scale
         check(H5LTset_attribute_double(gridGroup, ".", "Lx", &Lx, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Ly", &Ly, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Lz", &Lz, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Lv", &Lv, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Lm", &Lm, 1), DMESG("Attribute"));
         
         // Grid point number
         check(H5LTset_attribute_int(gridGroup, ".", "Nx", &Nx, 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nky", &Nky, 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nz", &Nz, 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nv", &Nv, 1), DMESG("Attribite"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nm", &Nm, 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Ns", &Ns, 1), DMESG("Attribute"));
         

 	// set Lengths
	check(H5LTset_attribute_double(gridGroup, ".", "X", &X(NxGlD), Nx), DMESG("Attribute"));
//	check(H5LTset_attribute_double(gridGroup, ".", "Y", &Y(NyGlD), Ny), DMESG("Attribute"));
	check(H5LTset_attribute_double(gridGroup, ".", "Z", &Z(NzGlD), Nz), DMESG("Attribute"));
	check(H5LTset_attribute_double(gridGroup, ".", "V", &V(NvGlD), Nv), DMESG("Attribute"));
	check(H5LTset_attribute_double(gridGroup, ".", "M", &M(NmGlD), Nm), DMESG("Attribute"));
         
    H5Gclose(gridGroup);

          
    }

     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

     int getGlobalSize() const { return Nx * Nky * Nz * Nv * Nm * Ns; };
     int getLocalSize() const { return NxLD * NkyLD * NzLD * NvLD * NmLD * NsLD;};
    
};



#endif // __GRID_H_

