#ifndef VLASOV_LORENTZ_H
#define VLASOV_LORENTZ_H

#include "Vlasov.h"
#include "FileAttr.h"
#include "Collisions.h"
	/* Helper function for collisional operator as defined in 
	*  PhD Thesis of Merz
	*/
class Vlasov_Lorentz : public Vlasov {
    
    /*
 Collisional frequency 
*/
protected:   
    
    double nu;

    

   
  public:


       Vlasov_Lorentz(Grid *grid, Parallel *parallel, Setup *setup, Geometry<HELIOS_GEOMETRY> *geo) : Vlasov(grid, parallel, setup, geo) {

        nu    = setup->get("Collision.Nu", 0.e0);
        
       };
 

       /**
        *   Set Data output parameters
        *
        * */
       void initDataOutput(hid_t fileID) {
     hid_t collisionGroup = check(H5Gcreate(fileID, "/Collisions",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Collision : H5Gcreate"));
     
     check(H5LTset_attribute_double(collisionGroup, ".", "Nu", &nu, 1), DMESG("H5LTset_attribute"));
     
     H5Gclose(collisionGroup);

  }   

       inline double Coll(const int x, const int y, const int z, const int v, const int m, const int s) {
          
          if(nu == 0.) return 0.;
          else {

              const double df1_dv    = (8. *(f1(x, y, z, v+1, m, s) - f1(x, y, z, v-1, m, s))  -1. *(f1(x, y, z, v+2, m, s) - f1(x, y, z, v-2, m, s)))/(12.*dv);
              const double ddf1_dvv  = (nu > 0.) ? (16. *(f1(x, y, z, v+1, m, s) + f1(x, y, z, v-1, m, s))  -1. *(f1(x, y, z, v+2, m, s) + f1(x, y, z, v-2, m, s)) - 30.*f1(x,y,z,v,m,s))/(12.*pow2(dv)): 0.;

              return    ((nu > 0.) ? 0.5 * nu * (ddf1_dvv - 2. * V(v) * df1_dv - pow2(V(v)) * ddf1_dvv): 0.);
          }
            return 0.;

       };
protected:
        /**
         * Program output
         *
         * */
        virtual void printOn(ostream &output) const {
         output   << "Collisions |  Lorentz " << ((nu > 0.) ? Num2String(nu *2.) :  "off") << std::endl;
        }
};




#endif // VLASOV_LORENTZ_H
