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
#include "external/SFL/SFL.h"

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
   __attribute__((vector)) static double sign(const double T) { return ((T >= 0.) ? 1. : -1); }
   /**
   * C has remainder modulo operator, we need real one
   **/ 
   int realmod(double x,double y) { return static_cast<int>(fmod((fmod(x,y) + y), y)); };

   /**
   *  @brief Bessel function of zeroth order
   *  @image html Bessel_J0.png
   *
   *  @depraciated
   *
   *  replaced by SFL
   **/
   //__attribute__((vector)) static inline double BesselJ0(const double x) { return j0(x); };

   /**
   *  @brief Modified bessel function of zeroth order
   *  @image html Bessel_I0.png
   *
   *  @todo any better approach ?
   *  
   *  @replaced by SFL
   */

   /**
   *  @brief Modified bessel function of first order
   *  @image html Bessel_I1.png
   *
   *  @todo any better approach ? Plenty of room for optimization
   *  @replaced by SFL
   **/


   /**
   *  @brief Modified bessel function of first order
   *  @image html Bessel_I1.png
   *
   **/
   static inline double _1mGamma0(const double b, bool gyro=true) {

    return gyro ? 1. - SFL::i0e(b) : b;

   };


   /**
   *  @brief Modified bessel function of first order
   *  @image html GK_1mG0_Pade.png
   *
   *  Defined as 
   *    \f[
   *        1 - \Gamma_{0,P} = b / \left( 1+ b \right).
   *    \f]
   *
   *  Note :
   *
   *  The corresponding drift-kinetic analog for $b \ll 1$ would be 
   *  $ 1 - \Gamma_{0,DK} = b $.
   *
   **/
   __attribute__((vector)) static inline double _1mGamma0_Pade(const double b) 
   {

    return  b / ( 1.e0 + b);

   };


   /**
   *  @brief Modified bessel function of first order
   *  @image html GK_Gamma0.png
   *
   *  Calculating I_0(b) exp(-b) individually, will overflow for
   *  I_0(b) due to exponetial behaviour for large b and underflow
   *  for exp(-b).
   *
   **/
   static inline double Gamma0(const double b, bool gyro=true) 
   {
     return gyro ? SFL::i0e(b) : b;
   };

   /**
   *  @brief Modified bessel function of first order
   *  @image html GK_Gamma1.png
   *
   *  Is this correct ?! Inconsistent with Asytrophys gyrokinetics
   *
   */
   static inline double Gamma1(const double b, bool gyro=true) 
   {
      return gyro ? SFL::i1e(b) : b;
   };

   /**
   *  @brief Modified bessel function of first order
   *  @image html G0mG1.png
   *
   */
   static inline double G0mG1(const double b, bool gyro=true) 
   {
      return gyro ? Gamma0(b) - Gamma1(b) : b;
   };

   /**
   *  @brief Modified bessel function of first order
   *  @image html Bessel_I1.png
   *  
   *  @note wtf ? 
   **/
   static inline double I1sp(const double b) 
   {
      //const double i = 1.;
      //return  -2. * i *  BesselI1(i * b)/b ;
      return - 2. * SFL::i1(b)/b;
   };



   /**
   *  @brief Modified bessel function of first order
   *   \f[     \Delta = \left( I_0(x) - I_1(x) \right) \exp(-x)  \f]
   *  @image html GK_Delta.png
   *
   **/
   static inline double Delta(const double x) 
   {
     return SFL::i0e(x) - SFL::i1e(x);

   }

   /**
   *  @brief Definition of \f$ \Delta_1 \f$  
   *   \f[   
   *       \int_0^\infty J_0^2(lambda) \mu  e^{-mu} d\mu = 
   *       \Delta(b) =  I_0(b) e^{-b} + b \left( I_1(b) - I_0(b) \right) e^{-b} 
   *    \f]
   *  @image html GK_Delta_1.png
   *
   **/
   static inline double Delta_1(const double b) 
   {
     return b * SFL::i1e(b) + (1.-b) * SFL::i0e(b);

   }







};

#endif // __SPECIAL_MATH_H
