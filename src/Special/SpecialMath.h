/*
 * =====================================================================================
 *
 *       Filename: Math.h
 *
 *    Description: Special math functions
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __SPECIAL_MATH_H
#define __SPECIAL_MATH_H

#include "Global.h"

/**
*  @brief Various math functions
*
*  @todo cleanup
*  @todo find more efficient ones
*
**/
class SpecialMath
{
  public:
   /**
   *   @brief returns sign of double
   *
   *
   **/ 
   static double sign(const double T) { return ((T >= 0.) ? 1. : -1); }
   /**
   * C has remainder modulo operator, we need real one
   **/ 
   int realmod(double x,double y) { return static_cast<int>(fmod((fmod(x,y) + y), y)); };

   /**
   *  @brief Bessel function of zeroth order
   *  @image html Bessel_J0.png
   *
   **/
   static inline double BesselJ0(const double x) { return j0(x); };

   /**
   *  @brief Modified bessel function of zeroth order
   *  @image html Bessel_I0.png
   *
   *  @todo any better approach ?
   *
   */
   static inline double BesselI0( const double x )
   {

     double ax,ans;

     if ((ax=fabs(x)) < 3.75) {
         const double y=pow2(x/3.75);
         ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
     } else {
         const double  y=3.75/ax;
         ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
          +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
          +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
          +y*0.392377e-2))))))));
     }
     
     return ans;
   }

   /**
   *  @brief Modified bessel function of first order
   *  @image html Bessel_I1.png
   *
   *  @todo any better approach ?
   **/
   static inline double BesselI1( double x ) 
   {

      double I1;
                                                           
      
      if (abs(x) <  3.75)  { 

        const double  p[8] = { 0.0, 0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3 };
        const double y=pow2((x/3.75));
        I1=x*(p[1]+y*(p[2]+y*(p[3]+y*(p[4]+y*(p[5]+y*(p[6]+y*p[7]))))));
      
      } else {
        
        const double  q[10] = { 0.0, 0.39894228e0, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 
                                         0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2 } ;
        const double ax=abs(x);
        const double y=3.75/ax;
        
        I1=(exp(ax)/sqrt(ax))*(q[1]+y*(q[2]+y*(q[3]+y*(q[4]+y*(q[5]+y*(q[6]+y*(q[7]+y*(q[8]+y*q[9]))))))));
        
      }
      
      return (x  > 0.) ? I1 : -I1;
   }



   /**
   *  @brief Modified bessel function of first order
   *  @image html Bessel_I1.png
   *
   **/
   static inline  double _1mGamma0(const double b, bool gyro=true) {

    if(gyro) return  1. - BesselI0(b)*exp(-b);
    else     return  b;

   };


   /**
   *  @brief Modified bessel function of first order
   *  @image html GK_1mG0_Pade.png
   *
   **/
   static inline  double _1mGamma0_Pade(const double b, bool gyro=true) 
   {

    if(gyro) return  b / ( 1.e0 + b);
    else     return  b;

   };


   /**
   *  @brief Modified bessel function of first order
   *  @image html GK_Gamma0.png
   *
   **/
   static inline  double Gamma0(const double b, bool gyro=true) 
   {

    if(gyro) return  BesselI0(b)*exp(-b);
    else     return  b;

   };

   /**
   *  @brief Modified bessel function of first order
   *  @image html GK_Gamma1.png
   *
   */
   static inline  double Gamma1(const double b, bool gyro=true) 
   {

      if(gyro) return  BesselI1(b)*exp(-b);
      else     return  b;

   };

   /**
   *  @brief Modified bessel function of first order
   *  @image html G0mG1.png
   *
   */
   static inline  double G0mG1(const double b, bool gyro=true) 
   {

      if(gyro) return  Gamma0(b) - Gamma1(b);
      else     return  b;

   };

   /**
   *  @brief Modified bessel function of first order
   *  @image html Bessel_I1.png
   *
   **/
   static inline  double I1sp(const double b) 
   {
	   const double i = 1.;
	   return  -2. * i *  BesselI1(i * b)/b ;
   };



   /**
   *  @brief Modified bessel function of first order
   *   \f[     \Delta = \left( I_0(x) - I_1(x) \right) \exp(-x)  \f]
   *  @image html GK_Delta.png
   *
   **/
   static inline  double Delta(const double x) 
   {

    return (BesselI0(x)-BesselI1(x))*exp(-x);

   }


};

#endif // __SPECIAL_MATH_H
