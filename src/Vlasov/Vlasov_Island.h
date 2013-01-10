
/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Aux.h
 *
 *    Description: Implementation of GK Vlasov's equation for two dimensional
 *                 slab geometries
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Global.h"

#ifndef __VLASOV_ISLAND_H
#define __VLASOV_ISLAND_H


#include "Vlasov_Aux.h"


/**
*  @brief Implementation of Vlasov's equation using Cilk Array Notation
*
*  Making extensively use of Intel's Cilk Plus Array Notation to faciliate
*  array operations (especially vectorization). (http://software.intel.com/en-us/articles/intel-cilk-plus/)
*  Supported by Intel(12.1)  and GCC (svn side branch)
*
*
**/
class VlasovIsland : public VlasovAux {

   nct::allocate ArrayX,
                 ArrayY;
  
   nct::allocate ArrayAp_mod; 

   CComplex *Ap_mod;
  
   CComplex  *G_lin,  ///< Modified phase distribution with fields contributions
            *Xi_lin; ///< Combined fields \f$ \Xi = \phi + v_\parallel A_{1\parallel}  \f$

   double width, ///< Magnetic Island width
          shear, ///< Shearing rate
          omega; ///< Island frequeny

   double *MagIs, *dMagIs_dx, *ky_filter;

   double p[3]; ///<  Coefficients for the island structure
 
   int  i   ; ///< Mode number of Island
   /**
   *    Please Document Me !
   *
   **/
   void  Vlasov_2D_Island(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinear                  [Nky][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs[NxGB], 
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const CComplex Ap_mod                    [NzLB][Nky][NxLB+4]      ,
                           CComplex Field0[Nq][NzLD][Nky][NxLD]   ,
                           const double dt, const int rk_step, const double rk[3]);

   void  Vlasov_2D_Island_EM(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex    nonLinearTerm               [Nky][NxLD  ][NvLD],
                           CComplex Xi       [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G        [NzLB][Nky][NxLB][NvLB],
                           CComplex Xi_lin       [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G_lin        [NzLB][Nky][NxLB][NvLB],
                           const CComplex Ap_mod                    [NzLB][Nky][NxLB+4]      ,
                           CComplex Field0[Nq][NzLD][Nky][NxLD]   ,
                           const double dt, const int rk_step, const double rk[3]);
   
   virtual void setupXiAndG_lin(
                           const CComplex g          [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4],
                                 CComplex Xi                     [NzLB][Nky][NxLB+4][NvLB],
                                 CComplex G                      [NzLB][Nky][NxLB  ][NvLB],
                           const int m, const int s);

   
   
   void  Vlasov_2D_Island_Rotation(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinear                  [Nky][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs[NxGB], 
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);
   
   void  Vlasov_2D_Island_Equi(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinear                  [Nky][NxLD  ][NvLD],
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);
   
   void  Vlasov_2D_Island_filter(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinear                  [Nky][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs[NxGB], 
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);


  public:

   /**
   *    Please Document Me !
   *
   **/
   VlasovIsland(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO * fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *bench, Collisions *coll); 
        
   /**
   *    Please Document Me !
   *
   **/
   void solve(std::string equation_tyoe, Fields *fields, CComplex  *fs, CComplex *fss, double dt, int rk_step, const double rk[3]);
 
  protected :
 
   void printOn(std::ostream &output) const;
   void initDataOutput(Setup *setup, FileIO *fileIO); 
 

};


#endif //  __VLASOV_ISLAND_H


