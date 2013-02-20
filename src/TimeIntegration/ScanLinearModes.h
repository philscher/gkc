/*
 * =====================================================================================
 *
 *       Filename: TimeIntegration.h
 *
 *    Description:Time Integration Interface for gkc. 
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#include "config.h"
#include "Global.h"

#ifndef SCANMODES_H__
#define SCANMODES_H__

#include "Setup.h"
#include "Grid.h"
#include "Parallel/Parallel.h"
#include "TimeIntegration.h"
#include "Fields.h"
#include "Init.h"
/**
*    @brief Class which handles explicit time intergration
*
*
*
**/
class ScanLinearModes  : public IfaceGKC {
 protected:
   
   Parallel *parallel;
   Fields   *fields;
   Vlasov   *vlasov;

   class ErrorVals
   {
      double *vals;
      int     size;
     public:

      ErrorVals(int _size) : size(_size) {
        vals = new double[size];
        for(int n = 0; n < size; n++) vals[0] = 0.;
      };
     ~ErrorVals() { delete[] vals; };

      void addValue(double val) {

         // push all elements down
         for(int n = 1; n < size-1; n++) vals[n] = vals[n-1];
         vals[0] = val;

      };
      
      double getMean()
      {
        double sum = 0.;
        for(int n = 0; n < size-1; n++) sum += vals[n];
        return sum/size;
      };

      /**
      *   @brief Calculate standard defiation
      *
      *   \f[
      *        \sigma = \sqrt{ \frac{1}{N} \sum_{i=1}^N left( x_i - \mu right)^2 }
      *   \f]
      *
      *
      **/
      double getStdDev()
      {
         double s2N = 0., mean = getMean();
         for(int n = 0; n < size-1; n++) s2N  += pow2(vals[n] - mean);

         return sqrt(s2N / size);
      };
   };

   std::vector<double> modes; ///< Stores the poloidal modes to be scanned
  
   double minAbsError_re, ///< Minumum absolute error (real frequency)
          minAbsError_im; ///< Minumum absolute error (growthrates)

   int InnerIterMax;   ///< Number of iterations before frequency is recalculated

   nct::allocate ArrayField0; ///< Array allocator for Field0 variable
   CComplex *Field0_t;
  
   /**
   *  @brief structure to store the eigenvalue
   *
   **/
   struct Frequency {
    double  ky       ;
    Complex frequency; ///< the eigenvalue
    Complex error    ; ///< the absolute error

   };

  TableAttr *freqTable;   ///< Eigenvalue Table
  hid_t     scanGroupID;  ///< HDF-5 reference to Eigenvalue group

 public:
   /**
   *   Constructor
   *
   **/ 
   ScanLinearModes(Setup *setup, Grid *grid, Parallel *parallel, Vlasov *vlasov, Fields *fields, FileIO *fileIO);
  ~ScanLinearModes();

   void solve(Vlasov *vlasov, Fields *fields, TimeIntegration *timeIntegration, Eigenvalue *eigenvalue, Init *init, Visualization *visual); 
   
   virtual void printOn(std::ostream &output) const {} ;
  
   void initData(Setup *setup, FileIO *fileIO);
};

#endif // SCANMODES_H__
