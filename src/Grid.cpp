/*
 * =====================================================================================
 *
 *       Filename: Grid.cpp
 *
 *    Description: Grid definitions including boundary
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Grid.h"
#include "Special/Integrate/Integrate.h"

/////////////////  Define Global Variables ////////////////////////
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
int    Nx, Nky, Nz, Nv, Nm, Ns, Nq;
double dx, dy, dz, dv,dt;

///////////////////////////////////////////////

Grid:: Grid (Setup *setup, Parallel *parallel, FileIO *fileIO) 
{
  
  muIntegrationType = setup->get("Grid.MuIntegrationType", "Gauss-Legendre");
    
  // Set Initial conditions
 
  // Global Domain Grid Number
  Nx  = setup->get("Grid.Nx", 32);
  Nky = setup->get("Grid.Nky", 9);
  Nz  = setup->get("Grid.Nz", 8 );
  Nv  = setup->get("Grid.Nv", 16);
  Nm  = setup->get("Grid.Nm", 1 );
  Ns  = setup->get("Grid.Ns", 1 );
  
  // Global Domain Grid Lengths
  Lx  = setup->get("Grid.Lx", 16.   );
  Ly  = setup->get("Grid.Ly", 32.   );
  Lz  = setup->get("Grid.Lz", 2. * M_PI);
  Lv  = setup->get("Grid.Lv", 4.  );
  Lm  = setup->get("Grid.Lm", 7.  );
  
  const double x0  = setup->get("Grid.x0", Lx/2.);
  //const double z0  = setup->get("Grid.z0", 0.   );
    
  int Ny = 2 * Nky - 2;
      
  check(((Nm == 1) && (Lm != 1.)) ? -1 : 0, DMESG("For Nm=1, set Lm = 1 !"));
  
  // Calculate grid distances values (equidistant for X, Y, Z, V)
  dx = (Nx > 1) ? Lx/((double) (Nx-1)) : Lx;
  dy = (Ny > 1) ? Ly/((double) (Ny-1)) : Ly;
  dz = (Nz > 1) ? Lz/((double) (Nz  )) : Lz;
  dv = (Nv > 1) ? 2.* Lv/ ( (double) (Nv-1)) : Lv;
      
  // Set Global Global lower/upper Domain/Boundary index
  NxGlD = 3; NxGuD=Nx+2; NxGlB=1; NxGuB=Nx+4; 
  NyGlD = 3; NyGuD=Ny+2; NyGlB=1; NyGuB=Ny+4; 
  NzGlD = 3; NzGuD=Nz+2; NzGlB=1; NzGuB=Nz+4; 
  NvGlD = 3; NvGuD=Nv+2; NvGlB=1; NvGuB=Nv+4; 
  NmGlD = 1; NmGuD=Nm;   NmGlB=1; NmGuB=Nm; 
  NsGlD = 1; NsGuD=Ns;   NsGlB=1; NsGuB=Ns; 
       
  NkyGlD = 0  ; NkyGuD=Nky-1; NkyGlB=0; NkyGuB=Nky-1; 
  NkyLlD = 0  ; NkyLuD=Nky-1; NkyLlB=0; NkyLuB=Nky-1;
  NkyLD  = Nky; NkyGD = Nky; 
  
  // Number of Ghost/Halo cells
  NxGC = 2; NyGC = 2; NzGC = 2; NvGC = 2; 

  //////////////////   Grid decomposition  ////////////////////////
  // Set local (L) lower (l) and upper (u) bounds
  // Currently we only have equal partition decomposition

  // Set local decomposition indices for X
  NxLlD = (NxGuD - NxGlD + 1)/parallel->decomposition[DIR_X] * parallel->Coord[DIR_X] + NxGC + 1;
  NxLuD = (NxGuD - NxGlD + 1)/parallel->decomposition[DIR_X] * (parallel->Coord[DIR_X]+1) + NxGC;
  NxLlB = NxLlD - NxGC; NxLuB = NxLuD + NxGC;

  // Set local decomposition indices for Y (wtf?)
  NyLlD = (NyGuD - NyGlD + 1)/parallel->decomposition[DIR_Y] *  parallel->Coord[DIR_Y] + NyGC + 1;
  NyLuD = (NyGuD - NyGlD + 1)/parallel->decomposition[DIR_Y] * (parallel->Coord[DIR_Y]+ 1) + NyGC;
  NyLlB = NyLlD - NyGC; NyLuB = NyLuD + NyGC;

  // Set local decomposition indices for Z
  NzLlD = (NzGuD - NzGlD + 1)/parallel->decomposition[DIR_Z] *  parallel->Coord[DIR_Z] + NzGC + 1;
  NzLuD = (NzGuD - NzGlD + 1)/parallel->decomposition[DIR_Z] * (parallel->Coord[DIR_Z]+1) + NzGC;
  NzLlB =  NzLlD - NzGC; NzLuB = NzLuD + NzGC;
    
  // Set local decomposition indices for V
  NvLlD = (NvGuD - NvGlD + 1)/parallel->decomposition[DIR_V] *  parallel->Coord[DIR_V] + NvGC + 1;
  NvLuD = (NvGuD - NvGlD + 1)/parallel->decomposition[DIR_V] * (parallel->Coord[DIR_V] + 1) + NvGC;
  NvLlB = NvLlD - NvGC; NvLuB = NvLuD + NvGC;
    
  // Set local decomposition indices for M (no boundary)
  NmLlD = (NmGuD - NmGlD + 1)/parallel->decomposition[DIR_M] *  parallel->Coord[DIR_M] + 1;
  NmLuD = (NmGuD - NmGlD + 1)/parallel->decomposition[DIR_M] * (parallel->Coord[DIR_M] + 1) ;
  NmLlB = NmLlD; NmLuB = NmLuD; 

  // Set local decomposition indices for S (no boundary)
  NsLlD = (NsGuD - NsGlD + 1)/parallel->decomposition[DIR_S] *  parallel->Coord[DIR_S] + 1;
  NsLuD = (NsGuD - NsGlD + 1)/parallel->decomposition[DIR_S] * (parallel->Coord[DIR_S] + 1) ;
  NsLlB = NsLlD; NsLuB = NsLuD; 

  // Set number of local domain points
  NxLD = NxLuD - NxLlD + 1; NyLD = NyLuD - NyLlD + 1;
  NzLD = NzLuD - NzLlD + 1; NvLD = NvLuD - NvLlD + 1;
  NmLD = NmLuD - NmLlD + 1; NsLD = NsLuD - NsLlD + 1;

  // Set number of local boundary points 
  NxLB = NxLuB - NxLlB + 1; NyLB = NyLuB - NyLlB + 1;
  NzLB = NzLuB - NzLlB + 1; NvLB = NvLuB - NvLlB + 1;
  NmLB = NmLuB - NmLlB + 1; NsLB = NsLuB - NsLlB + 1;
    
  // Set number of local/global Domain points 
  NxGB = Nx + 2*NxGC  ; NyGB = Ny + 2*NyGC;
  NzGB = Nz + 2*NxGC  ; NvGB = Nv + 2*NxGC;
  NmGB = Nm           ; NsGB = Ns;

  // Set Ranges of Local Domain (LD)
  RxLD.setRange(NxLlD,NxLD);
  RzLD.setRange(NzLlD,NzLD);
  RvLD.setRange(NvLlD,NvLD);
  RmLD.setRange(NmLlD,NmLD);
  RsLD.setRange(NsLlD,NsLD);
  RkyLD.setRange(NkyLlD,NkyLD);
    
  // Set Ranges of Local Boundary (LB)
  RxLB.setRange(NxLlB,NxLB);
  RzLB.setRange(NzLlB,NzLB);
  RvLB.setRange(NvLlB,NvLB);
  RmLB.setRange(NmLlB,NmLB);
  RsLB.setRange(NsLlB,NsLB);
    
  RxLB4.setRange(NxLlB-2,NxLB+4);
  
  // Set Ranges of Global Boundary (LB)
  RxGB.setRange(NxGlB, NxGB);
  RzGB.setRange(NzGlB, NzGB);
  RvGB.setRange(NvGlB, NvGB);
  RmGB.setRange(NmGlB, NmGB);
  RsGB.setRange(NsGlB, NsGB);
  
  // Set Ranges of Global Domain (LD)
  RxGD.setRange(NxGlD, Nx);
  RzGD.setRange(NzGlD, Nz);
  RvGD.setRange(NvGlD, Nv);
  RmGD.setRange(NmGlD, Nm);
  RsGD.setRange(NsGlD, Ns);
 
  ///////////////  Set Grid  Domain ////////////

  // X (Note : For gyro-averaged fields we have extended boundaries
  ArrayX = nct::allocate(nct::Range(NxGlB-2, NxGB+4))(&X);
  ArrayZ = nct::allocate(RzGB)(&Z);
  ArrayV = nct::allocate(RvGB)(&V);
  ArrayM = nct::allocate(RmGB)(&M, &dm);

  // Use  equidistant grid for X, Z and V
  bool includeX0Point = setup->get("Grid.IncludeX0Point", 0);
  for(int x = NxGlB-2; x <= NxGuB+2; x++) X[x] = -Lx/2. + dx * (x - NxGC - 1) + ((includeX0Point) ? dx/2. : 0.);
  //for(int z = NzGlB; z <= NzGuB; z++) Z[z] = -Lz/2. + dz * (z - NzGC - 1) ;
  for(int z = NzGlB; z <= NzGuB; z++) Z[z] =  dz * (z - NzGC - 1) ;
  for(int v = NvGlB; v <= NvGuB; v++) V[v] = -Lv + dv * (v - NvGlD);
    
  // M For mu we can choose between various integration type e.g. rectangle or Gaussian
  // @todo use setup to define integration method
  Integrate integrate(muIntegrationType, Nm, 0., Lm);

  for(int m = NmGlD, n = 0; m <= NmGuD; m++, n++) {

    M [m] = (Nm == 1) ? 0. : integrate.x(n) ; 
    dm[m] = (Nm == 1) ? 1. : integrate.w(n); 
  }

  // Set local Jacobian (better use geometry module ?)
  dXYZ  = dx * dy * dz;
  dXYZV = dx * dy * dz * dv;

  initData(fileIO);
}

Grid::~Grid () 
{

}

void Grid::printOn(std::ostream &output) const 
{

  bool doGyro = Nm > 1 ? true : false;
  
  output << "Domain     |  Lx : " << std::setw(4) << Lx << "   Ly : " << std::setw(4) << Ly 
         <<             "  Lz : " << std::setw(4) << Lz <<  "  Lv : " << std::setw(4) << Lv 
         << std::setw(4)          << (doGyro ? std::string( "  Lm : ") + Setup::num2str(Lm) : "") << std::endl
         << "Grid       |  Nx : " << std::setw(4) << Nx << "  Nky : " << std::setw(4) << Nky 
                     << "  Nz : " << std::setw(4) << Nz <<  "  Nv : " << std::setw(4) << Nv 
         << std::setw(4)          << (doGyro ? std::string( "  Nμ : ") + Setup::num2str(Nm) : "") << std::endl;
}
    
void Grid::initData(FileIO *fileIO) 
{

  hid_t gridGroup = fileIO->newGroup("/Grid");
  
  // set lengths
  check(H5LTset_attribute_double(gridGroup, ".", "Lx", &Lx, 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "Ly", &Ly, 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "Lz", &Lz, 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "Lv", &Lv, 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "Lm", &Lm, 1), DMESG("HDF-5 Error"));
         
  // set grid point number
  check(H5LTset_attribute_int(gridGroup, ".", "Nx", &Nx , 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_int(gridGroup, ".", "Nky",&Nky, 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_int(gridGroup, ".", "Nz", &Nz , 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_int(gridGroup, ".", "Nv", &Nv , 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_int(gridGroup, ".", "Nm", &Nm , 1), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_int(gridGroup, ".", "Ns", &Ns , 1), DMESG("HDF-5 Error"));
         
  // set grids 
  check(H5LTset_attribute_double(gridGroup, ".", "X", &X[NxGlD], Nx), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "Z", &Z[NzGlD], Nz), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "V", &V[NvGlD], Nv), DMESG("HDF-5 Error"));
  check(H5LTset_attribute_double(gridGroup, ".", "M", &M[NmGlD], Nm), DMESG("HDF-5 Error"));
  
  check(H5LTset_attribute_string(gridGroup, ".", "MuIntegrationType", muIntegrationType.c_str()), DMESG("HDF-5 Error")); 
  
  H5Gclose(gridGroup);
}

