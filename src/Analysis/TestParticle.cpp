
#include "TestParticle.h"


TestParticles::TestParticles(FileIO *fileIO, Setup *setup, Parallel *_parallel) : parallel(_parallel) 
{
    Total_Tracer = setup->get("Tracers.Number", 0);
    particles = new Particle[Total_Tracer];      
 


  /* 
  for(int n = 1; n <= NumberOfParticles; n++) {

         Particle->checkIsOnGrid;
         
         Particle->Advance();
         



      }

      // update Information   


 * */


   };


TestParticles::~TestParticles() {
//    delete FA_X;
//    delete FA_V;


}

void TestParticles::integrate(Vlasov *vlasov, Fields *fields, int step) {

   double x_pos = 0.;
   double y_pos = 0.;
   double z_pos = 0.;


   //note we need an interpolation scheme


     // Lorentz Equation \[ F = - q * (\nabla \phi - \left(\vec{v} \times\left( \vec{B} - \nabla A_{1\parallel} + B_{1\parallel} \right)  

   // interpolate derivative using tri-linear interpolation
   double dA1p_dx = 0., dA1p_dy = 0., dA1p_dz = 0., B1p = 0.;

   //dphi_dx = fields->phi(x,y,z,m,s);

   //F(x,y,z)  = p->charge * ( q + (velcocity * geo->B(x,y,z)))
      // update force
   double dphi_dx = 0.,  dphi_dy = 0.;

    for(int n = 0; n < Total_Tracer; n++) {
      /* 
        Particle p = particles[n];
      // check if position is on domain
       double dphi_dx, dphi_dy, dphi_dz;   
        const double F_x = - dphi_dx +    (p.v[DIR_Y] * B_z - p.v[DIR_Z] * B_y) + ( - p.v[DR_Z] * dA1p_dx                       ) + ( p.v[DIR_Y] * B1p);
      const double F_y = - dphi_dy +    (p.v[DIR_Z] * B_x - p.v[DIR_X] * B_z) + (   p.v[DIR_Z] * dA1p_dy                      ) + ( p.v[DIR_X] * B1p);
      const double F_z = - dphi_dz +    (p.v[DIR_Y] * B_z - p.v[DIR_Y] * B_x) + (   p.v[DIR_X] * dA1p_dx - p.v[DIR_Y] * dA1p_dy) + ( 0        );
       * */ 

      // update velocity
    }

}
    

void TestParticles::printOn(std::ostream &output) const {

             output   << "Tracer Particles    |  " << "Off" << std::endl;
        
}

  
void TestParticles::initData(Setup *setup, FileIO *fileIO) {

     hid_t particleGroup = check(H5Gcreate(fileIO->getFileID(), "/Tracer",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     
     hsize_t p_dim []      =  { Total_Tracer, 1};
     hsize_t p_mdim[]      =  { Total_Tracer, H5S_UNLIMITED };
     hsize_t p_offset[]    =  {  0 , 0};
     
     bool phiWrite = (parallel->Coord[DIR_VMS] == 0);
     
     
     FA_X      = new FileAttr("Postion" , particleGroup, 2, p_dim , p_mdim   , p_dim   , p_offset    ,  p_dim  , p_offset, true, fileIO->vector3D_tid);
     FA_V      = new FileAttr("Velocity" , particleGroup, 2, p_dim , p_mdim   , p_dim   , p_offset    ,  p_dim  , p_offset, true, fileIO->vector3D_tid);
     //FA_A      = new FileAttr("Position" , particleGroup, 4, p_dim , p_maxdim   , p_chunkdim   , p_moffset    ,  p_chunkBdim  , p_offset, true, fileIO->vector3D_tid);
     FA_Time  = fileIO->newTiming(particleGroup);


     H5Gclose(particleGroup);
     ////
     

  }   
     
  void TestParticles::writeData(Timing timing, double dt) 
{
    Vector3D X[Total_Tracer];
    Vector3D V[Total_Tracer];

    for(int n = 0; n < Total_Tracer; n++) {
      //X[n] = particles[n].p;
      //V[n] = particles[n].v;
    }
    
    FA_X->write(X);
    FA_V->write(V);
        
    parallel->print("Wrote Tracer data ... "); 

}
      
  void TestParticles::closeData() 
{
  
}
