/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Collision.h
 *
 *    Description: Handles signals and other triggers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef COLLISIONS__LORENTZ_H
#define COLLISIONS__LORENTZ_H

#include "Global.h"
#include "Setup.h"
#include "Geometry.h"
#include "Grid.h"

#include "FileAttr.h"

   /* Helper function for collisional operator as defined in 
   *  PhD Thesis of Merz
   */
class Collisions {
    
    /*
 Collisional frequency 
*/
public:   
    
    void calculate_BTerms(f) 
    {

        for(int m = NmLlD; m <=NmLuD; m++) {  for(int v = NvLlD; v <=NvLuD; v++) { 


            B1 += f1(x,y_k, z, v,m,s);
            B2 +=              f1(x,y_k, z, v,m,s);
            B3 +=              f1(x,y_k, z, v,m,s);
            B4 += pow2(V[v]) * f1(x,y_k, z, v,m,s);



        } }


    }


   // Pitch-Angle collisions
    void calculateCollisions() {

  
       Complex df_dm[NmLD];
       Complex df_dv[NmLD];

      // perform CD-4 derivative for dphi_dx , and dphi_dy
      for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k= NkyLlD; y_k <= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {
      for(int s=NsLlD; s<= NsLuD;s++) {
      
        for(int v = NvLlD; v <= NvLuD; v++) { for(int m  = NmLlD ;   m <= NmLuD ; m++  ) { 


                const Complex df_dv = (8.  *(f[s][m][z][y_k][x][v+1] - f[s][m][z][y_k][x][v-1]) - (f[s][m][z][y_k][x][v+2] - f[s][m][z][y_k][x][v-2]))/(12.*dv);
                df_dm[v][m] = (8.  *(f[s][m+1][z][y_k][x][v] - f[s][m-1][z][y_k][x][v]) - (f[s][m+2][z][y_k][x][v] - f[s][m-2][z][y_k][x][v]))/(12.*dm);


                C_v[v][m] =    (2. * plasma->B0 * M[m]) * df_dv -   ( 4. * M[m] * V[v]                 ) * df_dm 
                C_m[v][m] =  - (4. * M[m] * V[v]      ) * df_dv +   ( 8./plasma->B0 * pow2(V[v]) * M[m]) * df_dm
        } }


        for(int v = NvLlD; v <= NvLuD; v++) { for(int m  = NmLlD ;   m <= NmLuD ; m++  ) { 
        

                // use 4-th order central differences
                df_dt += ( 8. * (C_v[v+1][m] +  C_v[v-1][m]) - (C_v[v+2][m] - C_v[v-2][m])) / (12.*dv);
                df_dt += ( 8. * (C_m[v][m+1] +  C_v[v][m-1]) - (C_v[v][m+1] - C_v[v][m-1])) / (12.*dm);
      
        }}
      
      } } } }


    };


    double nu;

    Collisions(Setup *setup,  Grid *grid) {

        //vc = setup->get("Collisions.Frequency", 1.e-4);
        nu    = setup->get("Vlasov.Collision.Nu", 0.e0);
        double lambda_C = setup->get("Collisions.CoulombLogarithm", 1.e-4);
        double  e = 1.; 
        vc = M_PI * log(lambda_C) * pow4(e) * plasma->n_ref * plasma->L_ref / (cbrt(2.) * pow2(plasma->T_ref));


        pref_Cxy.resize(grid->RsGD); pref_Cxy = 0.;
        B.resize(grid->RsGD, Range(1,4)); B = 0.;
        // we need to calculate the pre-factors
        //

//        for(int s=NsGlD; s <= NsGuD; s++)  pref_Cxy(s) += - vc * plasma->species(j).n * pow2(plasma->species(j).q)* plasma->species(j).T * pow(plasma->species[s].m, 3./2.) 



        // caluclate B's for conservation term
     /* 
        double v2 = pow2(V[v]) + M[m]; 
 
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        :1
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);

      * */ 


        // normalized velocity : double x = v/scale_v;
    };

// will probably be moved to Vlasov soon, why not by multiderivative , so I dont't have to mess
//  up the interface and it is still clean, also need some boundaries ?
    /**
    * Derivative of the error function 
     * By definition
     *  \f[ Derf(x) = \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}} \exp\left(-x^2 \right)  \f]
     * */
   inline double Derf(const double x) { return 2./sqrt(M_PI)*exp(-pow2(x)); };


    /**
     *  Calculates helper functions
     *  \f[ F_1(x) = x \frac{d erf(x)}{dx} + \left( 2 x^2 - \right) erf(x) \f]
     * */
   inline double F1(const double x) { return x*Derf(x) + (2. * pow2(x) - 1.) * erf(x); };
    /**
     *  Calculates helper functions
     *  \f[ F_2(x) = \left( 1 - \frac{2}{3} x^2 \right) erf(x) - x \frac{d erf}{dx} \f]
     * */
   inline double F2(const double x) { return (1.-2./3.*pow2(x)) * erf(x) - x* Derf(x); };
    /**
     *  Calculates helper functions
     *  \f[ F_2(x) = F_1(x) + 3 F_2(x) \f]
     * */
   inline double F3(const double x) { return F1(x)+3.*F2(x); };

    
   

    /**
     *    \f[ \left< C_{\sigma j}^{\perp} \right> = - \frac{1}{v^5} \left( 2 v^2 F_1 + 3 B_0 \mu F_2  \right) \nabla_\perp f_\sigma \f]
     *
     **/
    inline double C_xy(const double df_dp, const int x, const int y, const int z, const int v, const int m, const int s) {
        
      const double xc = 0.;


       return pref_Cxy(s) * 1./pow5(V[v]) * ( V[v] * F1(xc) + 3. * plasma->B0 * M[m] * df_dp / plasma->species[s].m * F2(df_dp));

    };


    /**
     *    \f[ \left< C_{\sigma j}^{\perp} \right> = - \frac{1}{v^5} \left( 2 v^2 F_1 + 3 B_0 \mu F_2  \right) \nabla_\perp f_\sigma \f]
     *
     **/
    inline double C_vm(const double df_dv, const double df_dm, const double f, const int x, const int y, const int z, const int v, const int m, const int s) {
  
      int j = s;
       const double xc = plasma->species[s].scale_v / plasma->species(j).scale_v * V[v];
      
       // parallel velocity
      const double C_v =    ((2. * plasma->B0 * M[m])/ plasma->species[s].m * F1(xc) + pow2(V[v]) * F3(x)) * df_dv + (6. * M[m] * pow2(V[v]) * F2(xc)) * df_dm - 2. * F3(xc)/pow3(V[v]) * V[v] * f;
      const double C_m =    ((6. * M[m] * V[v] * F2(xc)) * df_dv  + 4./plasma->B0 * pow2(V[v]) * M[m] * F1(xc)) * df_dv + (4. * pow2(M[m]) * F3(xc)) * df_dm - 2. * F3(xc)/pow3(V[v]) * 2. * M[m] * f;
        
       
      // BUG : this is wrong, we need to diffrentiate over v
      const double Cvm = (C_v + C_m) / dv;

       return Cvm;

    };


   inline double A() { return 0.0;};

  public:
       Collisions(Setup *setup, Geometry *geometry) {};
  
       void initDataOutput(hid_t fileID) {
  //######################################################################################
     // For electric and magnetic potential
     hid_t collisionGroup = check(H5Gcreate(fileID, "/Collisions",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     
     check(H5LTset_attribute_double(collisionGroup, ".", "Nu", &nu, 1), DMESG("H5LTset_attribute"));
     
     H5Gclose(collisionGroup);

  }   
     
protected:

        virtual void printOn(std::ostream &output) const {
         output   << "Collisions |  " << ((nu > 0.) ? Num2String(nu) :  "off") << std::endl;
        }
};




#endif // COLLISIONS_H
