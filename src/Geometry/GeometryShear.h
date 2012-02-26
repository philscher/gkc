/*
 * =====================================================================================
 *
 *       Filename: GeometryShear.h
 *
 *    Description: Definition of shearless slab geomtry
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_SHEARED_H
#define GEOMETRY_SHEARED_H

#include "Geometry.h"
#include "Global.h"
#include "Setup.h"

#include "FileIO.h"

/**
 *
 *    Gives the Geometric coefficients for a sheared slab geometry with
 *    a magnetic field accroding to.
 *    \f[ \vec{B}_0 / B_0 = \vec{b} = \left(0,-x/L_s,1 \right) \f].
 *
 *    Reference : Jenko, PhD 2001, Section 3.5 Flussschlauchgeometry
 *
 *    where $L_s$ is defined as the connection length. The metric takes
 *    the simple form 
 *
 *    \f[
 *          g^{ij} = 
 *              \left( \begin{array}{lll}
 *                      1     &  z/L_s                 & 0 \\
 *                       z\L_s &  1 + \frac{z^2}{L_s^2} & 0 \\
 *                       0     &  0                     & 1
 *                       \end{array} \right)
 *    \f]
 *
 *    This results in the spatial operators of the form
 *
 *    \f[ \vec{b} \cdot \nabla = \partial_z \f]
 *    \f[ \nabla_\perp^2 = \frac{\partial^2}{\partial x^2} + \left(1 + \frac{z^2}{L_s^2} \right)
 *              \frac{\partial^2}{\partial y^2} + 2 \frac{z}{L_s} \frac{\partial^2}{\partial x \partial y} \f]
 *    \f[ \vec{b} \times \nabla A \cdot \nabla = \frac{\partial A}{\partial x}{\partial}{\partial y} -
 *                  \frac{\partial A}{\partial y}{\partial }{\partial x}  \f]
 *
 *  also, with
 *  L_s = 
 *  L_c = 2 \pi q R
 *  L_s =  \frac{q R}{\hat{s}}
 *
 *  with \hat{s} the shear.
 *  and for the parallel length L_z, we need to choose the connection length L_c. So once L_s and L_z
 *  is given we calculate the shear to
 *  \f[ \hat{s} = \frac{L_z }{2 \pi L_s} \f]
 *  
 *  for consistency  with the perdiocic parallel boundary conditions we ned to fullfill the 
 *  relation
 *  \f[ 2 \pi \hat{s} L_x = n_s L_y  or \frac{L_z}{L_s} = n * L_y \f]
 *  where n_s is an integer number, or an interpolation method needs to be used. 
 *
 *
 *  ToDO : Parallelize the parallel boundary condition.
 *
 *
 * */


class GeometryShear : public Geometry<GeometryShear>
{
  double Ls;
  double shear;

  bool connectFieldLines, roundShearToConnect;
 public:

   void printOn(ostream& output) const {
         output   << "Geometry  |  Sheared Slab   shear : " << shear << " (Ls : " << Ls << " ) "  << " eps_hat : " << eps_hat << std::endl
                  << "          |  Connect Field Lines (" << (connectFieldLines ? "true" : "false") << ")   RoundShear (" << (roundShearToConnect ? "true" : "false") << ")" << std::endl;
   };



  GeometryShear(Setup *setup, FileIO *fileIO) : Geometry<GeometryShear>(setup, fileIO,  true)  {
     

    shear               = setup->get("Geometry.Shear"   , 0.);
    connectFieldLines   = setup->get("Geometry.ConnectFieldLines", 1);
    roundShearToConnect = setup->get("Geometry.RoundShearToConnect", 1);

    //check connection length, underwise parallel boundary will fail
    // modf :  extract signed integral and fractional values from floating-point number
    {
        if((abs(fmod(2. * M_PI * shear * Lx/Ly, 1.)) > 1.e-5) && roundShearToConnect) {
		std::cout << "Warning : Rounding shear to ensure magnetic field line connection" << std::endl;
                const double di = round(2. * M_PI * shear * Lx/Ly);
  		shear = di / (2. * M_PI) * Ly/Lx;
	}

    }
    
    // if shear is zero, set connection length to high number (what about FPE)
    Ls = (shear != 0.) ?   Lz / (2. * M_PI * shear) : 1.e99;
  };

  
  inline double J(const int x, const int y, const int z) { return 1.; };



  // define metric elements
  inline double g_xx(const int x, const int y, const int z) { return 1.0;                 };
  inline double g_xy(const int x, const int y, const int z) { return Z(z)/Ls;             }; 
  inline double g_xz(const int x, const int y, const int z) { return 0.0;                 };
  inline double g_yy(const int x, const int y, const int z) { return 1.0 + pow2(Z(z)/Ls); };
  inline double g_yz(const int x, const int y, const int z) { return 0.0;                 };
  inline double g_zz(const int x, const int y, const int z) { return 1.0;                 };
  
  // define magnetic field and its variations
  inline double B      (const int x, const int y, const int z) { return 1.; };
  
  inline double dB_dx  (const int x, const int y, const int z) { return 0.; };
  inline double dB_dy  (const int x, const int y, const int z) { return 0.; };
  inline double dB_dz  (const int x, const int y, const int z) { return 0.; };
 
  inline double Kx(const int x, const int y, const int z)  { return 0.; };
  inline double Ky(const int x, const int y, const int z)  { return 0.; };
  
  /**
   *
   *   Get the value of shear at position x, which is constant for sheared slab geometry 
   *
   */
  inline double getShear(const int x)  { return shear; };
  


  // for sheared magnetic fields, we have special boundary conditions, because
  // we need to take care of possible shear
  inline ShearB getYPos(const int x,const int y) {

//            // calculte shift in y, 
//            // then map by real modulo (add ghost cells)
            const double sy = 2. * M_PI * shear * X(x);
            const int yi = ((int) (sy/dy));

         
            int lDy, uDy;
            if(connectFieldLines) {
//	        lDy = (NzLlD == NzGlD) ? realmod(y- 3 - yi, Ny) + 3: y;
//                uDy = (NzLuD == NzGuD) ? realmod(y-3 + yi, Ny) + 3: y;
	        
            } 
	    else {
		lDy = y; uDy = y;
            }
            // if field line leaves the domain, don't connect them anymore
       /*     else {
	            lDy = (NzLlD == NzGlD) ? ((y-3 + yi > Ny) ? -1 : y + yi + 3) : y;
                uDy = (NzLuD == NzGuD) ? ((y-3 - yi <  0) ? -1 : y - yi + 3) : y;
            }
         */   
            return ShearB(lDy, uDy, 0.0  );//(yv - Y(yi-1))/dy);
              
            check(-1, DMESG("Should not be reached here"));        
    };


    std::string getGeometryName() { return "Sheared Slab"; };

    void initDataOutput(FileIO *fileIO, hid_t geometryGroup) {
          check(H5LTset_attribute_string(geometryGroup, ".", "Type", "Sheared Slab"), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_string(geometryGroup, ".", "ConnectFieldLines", connectFieldLines ? "true" : "false"), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_string(geometryGroup, ".", "RoundShearToConnect", roundShearToConnect ? "true" : "false"), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "Shear"   ,  &shear, 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "eps_hat"   ,  &eps_hat, 1), DMESG("H5LTset_attribute"));
    }


// transformation
  double getY(const int x, const int y, const int z) {
	return Y(y)+X(x)*Z(z)/Ls;
 }
//   void initDataOutput(FileIO *fileIO) {
//   } ;
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

private:

  
};


#endif // GEOMETRY_H
