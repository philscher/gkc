/*
 * =====================================================================================
 *
 *       Filename: Grid.cpp
 *
 *    Description: Grid definitions incuding boundary
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Grid.h"

#include "Special/Integrate/Integrate.h"


// **************** Define Global Variables ************* //
int NxLlD, NxLuD, NxLlB, NxLuB; 
int NyLlD, NyLuD, NyLlB, NyLuB; 
int NkyLlD, NkyLuD, NkyLlB, NkyLuB; 
int NzLlD, NzLuD, NzLlB, NzLuB; 
int NvLlD, NvLuD, NvLlB, NvLuB; 
int NmLlD, NmLuD, NmLlB, NmLuB; 
int NsLlD, NsLuD, NsLlB, NsLuB; 
      
int NxLD, NyLD, NkyLD, NzLD, NvLD, NmLD, NsLD;
int NxLB, NyLB, NzLB, NvLB, NmLB, NsLB;
int NxGB, NyGB, NkyGD, NzGB, NvGB, NmGB, NsGB;

int KxLuD, KyLuD, KzLuD, KxLlD, KyLlD, KzLlD;

int NkyGlD, NkyGuD, NkyGlB, NkyGuB; 

int NxGlD, NxGuD, NxGlB, NxGuB; 
int NyGlD, NyGuD, NyGlB, NyGuB; 
int NzGlD, NzGuD, NzGlB, NzGuB; 
int NvGlD, NvGuD, NvGlB, NvGuB; 
int NmGlD, NmGuD, NmGlB, NmGuB; 
int NsGlD, NsGuD, NsGlB, NsGuB; 

double *X, *Y, *V, *M, *Z;

  double Lx, Ly, Lz, Lv, Lm;
  int    Nx, Nky, Nz, Nv, Nm, Ns;
  double dx, dy, dz, dv,dt;

// ************************************************* //





Grid:: Grid (Setup *setup, Parallel *parallel, FileIO *fileIO) 
{
      // Set Initial conditions
      Nx  = setup->get("Grid.Nx", 32);
      Nky = setup->get("Grid.Nky", 9);
      Nz  = setup->get("Grid.Nz", 8 );
      Nv  = setup->get("Grid.Nv", 16);
      Nm  = setup->get("Grid.Nm", 1 );
      Ns  = setup->get("Grid.Ns", 1 );
      Lx  = setup->get("Grid.Lx", 16.   );
      Ly  = setup->get("Grid.Ly", 32.   );
      Lz  = setup->get("Grid.Lz", 2. * M_PI);
      Lv  = setup->get("Grid.Lv", 4.  );
      Lm  = setup->get("Grid.Lm", 7.  );
    
      int Ny = 2 * Nky - 2;
      
      check(((Nm == 1) && (Lm != 1.)) ? -1 : 0, DMESG("For Nm=1, set Lm = 1 !"));

      // Calculate some preknown values
       dx = (Nx > 1) ? Lx/((double) (Nx-1)) : Lx;
       dy = (Ny > 1) ? Ly/((double) (Ny-1)) : Ly;
 //      dz = (Nz > 1) ? Lz/((double) (Nz-1)) : Lz;
       dz = (Nz > 1) ? Lz/((double) (Nz  )) : Lz;
       //dm = (Nm > 1) ? Lm/((double) (Nm-1)) : 1.;
       dv = (Nv > 1) ? 2.* Lv/ ( (double) (Nv-1)) : Lv;
      
      /* 
       dy = Ly     /  (double) Ny;
       dz = Lz     /  (double) Nz;
       dv = 2. * Lv/  (double) Nv;
       dm = (Nm > 1) ? Lm/ (double) Nm : 1.e0;
       */
       NxGlD = 3; NxGuD=Nx+2; NxGlB=1; NxGuB=Nx+4; 
       NyGlD = 3; NyGuD=Ny+2; NyGlB=1; NyGuB=Ny+4; 
       NzGlD = 3; NzGuD=Nz+2; NzGlB=1; NzGuB=Nz+4; 
       NvGlD = 3; NvGuD=Nv+2; NvGlB=1; NvGuB=Nv+4; 
       NmGlD = 1; NmGuD=Nm;   NmGlB=1; NmGuB=Nm; 
       NsGlD = 1; NsGuD=Ns;   NsGlB=1; NsGuB=Ns; 
       
       NkyGlD = 0; NkyGuD=Nky-1; NkyGlB=0; NkyGuB=Nky-1; 
       NkyLlD = 0; NkyLuD=Nky-1; NkyLlB=0; NkyLuB=Nky-1;
       NkyLD=Nky; NkyGD=Nky; NkyGD=Nky;

    // Set Domain for all processes

       // Number of Ghost/Halo cells
       NxGC = 2 ;    NyGC = 2 ;    NzGC = 2 ;    NvGC = 2 ; 

    // Set local (L) lower (l) and upper (u) bounds
    NxLlD = (NxGuD - NxGlD + 1)/parallel->decomposition[DIR_X] * parallel->Coord[DIR_X] + NxGC + 1;
    NxLuD = (NxGuD - NxGlD + 1)/parallel->decomposition[DIR_X] * (parallel->Coord[DIR_X]+1) + NxGC;
    NxLlB=NxLlD - NxGC; NxLuB=NxLuD + NxGC;

    NyLlD = (NyGuD - NyGlD + 1)/parallel->decomposition[DIR_Y] * parallel->Coord[DIR_Y] + NyGC + 1;
    NyLuD = (NyGuD - NyGlD + 1)/parallel->decomposition[DIR_Y] * (parallel->Coord[DIR_Y]+1) + NyGC;
    NyLlB=NyLlD - NyGC; NyLuB=NyLuD + NyGC;

    NzLlD = (NzGuD - NzGlD + 1)/parallel->decomposition[DIR_Z] * parallel->Coord[DIR_Z] + NzGC + 1;
    NzLuD = (NzGuD - NzGlD + 1)/parallel->decomposition[DIR_Z] * (parallel->Coord[DIR_Z]+1) + NzGC;
    NzLlB =  NzLlD - NzGC; NzLuB=NzLuD + NzGC;
    
    // Settings for V
    NvLlD = (NvGuD - NvGlD + 1)/parallel->decomposition[DIR_V] *  parallel->Coord[DIR_V] + NvGC + 1;
    NvLuD = (NvGuD - NvGlD + 1)/parallel->decomposition[DIR_V] * (parallel->Coord[DIR_V] + 1) + NvGC;
    NvLlB=NvLlD - NvGC; NvLuB=NvLuD + NvGC;
    
    // Settings for S (no boundary cells)
    NmLlD = (NmGuD - NmGlD + 1)/parallel->decomposition[DIR_M] *  parallel->Coord[DIR_M] + 1;
    NmLuD = (NmGuD - NmGlD + 1)/parallel->decomposition[DIR_M] * (parallel->Coord[DIR_M] + 1) ;
    NmLlB = NmLlD; NmLuB = NmLuD; 


    // Settings for S (no boundary cells)
    NsLlD = (NsGuD - NsGlD + 1)/parallel->decomposition[DIR_S] *  parallel->Coord[DIR_S] + 1;
    NsLuD = (NsGuD - NsGlD + 1)/parallel->decomposition[DIR_S] * (parallel->Coord[DIR_S] + 1) ;
    NsLlB = NsLlD; NsLuB = NsLuD; 


    
    NxLD = NxLuD - NxLlD + 1; NyLD = NyLuD - NyLlD + 1;
    NzLD = NzLuD - NzLlD + 1; NvLD = NvLuD - NvLlD + 1;
    NmLD = NmLuD - NmLlD + 1; NsLD = NsLuD - NsLlD + 1;

    NxLB = NxLuB - NxLlB + 1; NyLB = NyLuB - NyLlB + 1;
    NzLB = NzLuB - NzLlB + 1; NvLB = NvLuB - NvLlB + 1;
    NmLB = NmLuB - NmLlB + 1; NsLB = NsLuB - NsLlB + 1;
    
    NxGD = Nx ; NyGD = Ny; NzGD = Nz; NvGD = Nv; NmGD= Nm; NsGD=Ns;
    NxGB = Nx + 2*NxGC  ; NyGB = Ny + 2*NyGC;
    NzGB = Nz + 2*NxGC  ; NvGB = Nv + 2*NxGC;
    NmGB = Nm           ; NsGB = Ns;

    // Ranges of Local Domain (LD)
    RxLD.setRange(NxLlD,NxLD);
    //RyLD.setRange(NyLlD,NyLD);
    RzLD.setRange(NzLlD,NzLD);
    RvLD.setRange(NvLlD,NvLD);
    RmLD.setRange(NmLlD,NmLD);
    RsLD.setRange(NsLlD,NsLD);
    
    RkyLD.setRange(NkyLlD,NkyLD);
  
    // Ranges of Local Boundary (LB)
    RxLB.setRange(NxLlB,NxLB);
    //RyLB.setRange(NyLlB,NyLB);
    RzLB.setRange(NzLlB,NzLB);
    RvLB.setRange(NvLlB,NvLB);
    RmLB.setRange(NmLlB,NmLB);
    RsLB.setRange(NsLlB,NsLB);
    
    RxLB4.setRange(NxLlB-2,NxLB+4);
   // RyLB4.setRange(NyLlB-2,NyLB+4);
  
    // Ranges of Global Boundary (LB)
    RxGB.setRange(NxGlB, NxGB);
   // RyGB.setRange(NyGlB, NyGB);
    RzGB.setRange(NzGlB, NzGB);
    RvGB.setRange(NvGlB, NvGB);
    RmGB.setRange(NmGlB, NmGB);
    RsGB.setRange(NsGlB, NsGB);
  
    // Ranges of Global Domain (LD)
    RxGD.setRange(NxGlD, NxGD);
  //  RyGD.setRange(NyGlD, NyGD);
    RzGD.setRange(NzGlD, NzGD);
    RvGD.setRange(NvGlD, NvGD);
    RmGD.setRange(NmGlD, NmGD);
    RsGD.setRange(NsGlD, NsGD);
 
    // Resize Send/Recv buffers for program to fit the ghost cells
    //RB.setRange(0,2);
    //RB4.setRange(0,4);


    ///////////////  Set Grid  Domain ////////////

   // Init Grid variables (-2, 2) extra points
    //X.resize(Range(NxGlB-2, NxGuB+2));
    
    // X
    nct::allocate(nct::Range(NxGlB-2, NxGB+4))(&X);
    bool includeX0Point = setup->get("Grid.IncludeX0Point", 0);
    for(int x = NxGlB-2; x <= NxGuB+2; x++) X[x] = -  Lx/2. + dx * ( x - NxGC - 1) + ((includeX0Point) ? dx/2. : 0.);


    // Use  equidistant grid for Y, Z and V
    nct::allocate(nct::Range(NyGlB  , NyGB))(&Y);
    for(int y = NyGlB; y <= NyGuB; y++) Y[y] = dy * ( y - NyGC - 1);
    nct::allocate(nct::Range(NzGlB  , NzGB))(&Z);
    for(int z = NzGlB; z <= NzGuB; z++) Z[z] = dz * ( z - NzGC - 1);
    nct::allocate(nct::Range(NvGlB  , NvGB))(&V);
    for(int v = NvGlB; v <= NvGuB; v++) V[v] = -  Lv + dv * ( v - NvGlD);
    
    // M For mu we can choose between linear and Gaussian integration
    // In case of Nm=1 (drift-kinetic, gyro-1st), dm = 1, and M = 0. ! shoudl be automatic ?
    nct::allocate(nct::Range(NmGlB  , NmGB))(&M, &dm);
    Integrate integrate("Gauss-Legendre", Nm, 0., Lm);
    for(int m=NmGlD, n=0; m <= NmGuD; m++, n++) { M[m] = (Nm == 1) ? 0. : integrate.x(n) ; dm[m] = (Nm ==1) ? 1. : integrate.w(n); }

    //
   dXYZ  = dx * dy * dz;
   dXYZV = dx * dy * dz * dv;

    initData(fileIO);

}

Grid::~Grid () {

}




void Grid::printOn(std::ostream &output) const {
    bool do_gyro = Nm > 1 ? true : false;
         output   << "Domain     |  Lx : " << Lx << "  Ly : " << Ly << "  Lz : " << Lz << "  Lv : " << Lv << ((do_gyro) ? std::string("  Lm : ") + Setup::num2str(Lm) : "") << std::endl
                  << "Grid       |  Nx : " << Nx << "  Nky : " << Nky << "  Nz : " << Nz << "  Nv : " << Nv << ((do_gyro) ? std::string("  Nμ : ") + Setup::num2str(Nm) : "") << std::endl;
     
            
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
    
void Grid::initData(FileIO *fileIO) {
          
      hid_t gridGroup = check(H5Gcreate(fileIO->getFileID(), "/Grid",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
          

         
          // Length scale
         check(H5LTset_attribute_double(gridGroup, ".", "Lx", &Lx, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Ly", &Ly, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Lz", &Lz, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Lv", &Lv, 1), DMESG("Attribute"));
         check(H5LTset_attribute_double(gridGroup, ".", "Lm", &Lm, 1), DMESG("Attribute"));
         
         // Grid point number
         check(H5LTset_attribute_int(gridGroup, ".", "Nx", &Nx , 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nky",&Nky, 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nz", &Nz , 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nv", &Nv , 1), DMESG("Attribite"));
         check(H5LTset_attribute_int(gridGroup, ".", "Nm", &Nm , 1), DMESG("Attribute"));
         check(H5LTset_attribute_int(gridGroup, ".", "Ns", &Ns , 1), DMESG("Attribute"));
         

    // set Lengths
   check(H5LTset_attribute_double(gridGroup, ".", "X", &X[NxGlD], Nx), DMESG("Attribute"));
//   check(H5LTset_attribute_double(gridGroup, ".", "Y", &Y(NyGlD), Ny), DMESG("Attribute"));
   check(H5LTset_attribute_double(gridGroup, ".", "Z", &Z[NzGlD], Nz), DMESG("Attribute"));
   check(H5LTset_attribute_double(gridGroup, ".", "V", &V[NvGlD], Nv), DMESG("Attribute"));
   check(H5LTset_attribute_double(gridGroup, ".", "M", &M[NmGlD], Nm), DMESG("Attribute"));
         
    H5Gclose(gridGroup);

          
    }
