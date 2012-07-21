#ifndef VLASOV_LENARD_BERNSTEIN_H
#define VLASOV_LENARD_BERNSTEIN_H

#include "Vlasov.h"
#include "Collisions.h"



/**
 *  \defgroup Vlasov_LenardBernstein '''Lenard-Bernstein collision model''
 *
 *  \sa 
 *     A. Lenard and Ira B. Bernstein, Plasma Oscillations with Diffusion in Velocity space,
 *     Physical Review Letters (1958)
 *     
 *   Small angel collision with
 *
 *   \f[ \mathcal{C}{()} =  \beta \frac{\partial}{\partial v} \left(  v F_{1}
 *                                  + v^2_0 \frac{\partial F_1}{\partial v} \right) \f]
 *   with \f[ v_0 \f], the root mean square speed corresponding to the equilibrium distribution \f[ F_0 \f]
 *  Quare of Root mean square defined as \f[ v_{rms}^2 = \int_0^{\infty} v^2 F_{0\sigma}  = 1/2\f] and
 *  Folowing Lenard and Bernstein, \f[ \beta \f] can be roughly estimated by comapring with the true
 *  Fokker-Planck equation. One get's approximately \f[ \beta \approx \frac{4 \ pi e^4}{m^2 v_0^3} \f].
 *
 *  \
 */
class Vlasov_LenardBernstein : public Vlasov {
    
    /*
 Collisional frequency 
*/
protected:   
    
    double beta;

   
  public:


       Vlasov_LenardBernstein(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry<GKC_GEOMETRY> *geo, FFTSolver *fft) : Vlasov(grid, parallel, setup, fileIO, geo, fft) {

        beta    = setup->get("Collision.Beta", 0.e0);
        
       };
 

       /**
        *   Set Data output parameters
        *
        * */
       virtual void initDataOutput(hid_t fileID) {
            hid_t collisionGroup = check(H5Gcreate(fileID, "/Collisions",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Collision : H5Gcreate"));
     
            check(H5LTset_attribute_string(collisionGroup, ".", "Model", "Lenard-Bernstein"), DMESG("H5LTset_attribute"));
            check(H5LTset_attribute_double(collisionGroup, ".", "Beta", &beta, 1), DMESG("H5LTset_attribute"));
            H5Gclose(collisionGroup);

  }   

       inline double Coll(const int x, const int y, const int z, const int v, const int m, const int s) {
          
          if(beta == 0.) return 0.;
          else {

              //const double v2_rms = 3./2. * pow2(plasma->species(s).scale_v);
              const cmplxd v2_rms = pow2(plasma->species(s).scale_v);

  //            const cmplxd df1_dv    = (8. *(f1(x, y, z, v+1, m, s) - f1(x, y, z, v-1, m, s))  -1. *(f1(x, y, z, v+2, m, s) - f1(x, y, z, v-2, m, s)))/(12.*dv);
//              const cmplxd ddf1_dvv  = (beta > 0.) ? (16. *(f1(x, y, z, v+1, m, s) + f1(x, y, z, v-1, m, s))  -1. *(f1(x, y, z, v+2, m, s) + f1(x, y, z, v-2, m, s)) - 30.*f1(x,y,z,v,m,s))/(12.*pow2(dv)): 0.;

          //    return    ((beta > 0.) ? beta * (f1(x,y,z,v,m,s)  + V(v) * df1_dv + v2_rms * ddf1_dvv): 0.);
          }
            return 0.;

       };
protected:
        /**
         * Program output
         *
         * */
        virtual void printOn(ostream &output) const {
         output   << "Collisions |  Model : Lenard-Bernstein    Beta = " << beta << std::endl;
        }
};




#endif // VLASOV_LENARD_BERNSTEIN_H
