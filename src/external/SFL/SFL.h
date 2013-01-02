
#ifndef FFN__H_
#define FFN__H_


/**  
*    The Special Function Library.
*
*   
*   Original File from :  http:  *people.sc.fsu.edu/~jburkardt/cpp_src/fn/fn.html
*
*   included following notes
*
*
*   (c) This code is distributed under the GNU LGPL license version 3 or 
*       at your opinion any later version
*
*   The class structure is distributed under the BSD license. However, each
*   function has an individual license, please refer to it.
*
*  @note  FN is a C++ library which contains a selection of routines for the evaluation of 
*   certain elementary and special functions, by Wayne Fullerton.
*
*   The original version of the library provided routines for complex, single precision real,
*   and double precision real arguments and used the prefixes "C" and "D" to indicate the complex 
*   and double precision versions.
*
*   This scheme has been modified for consistency, and also to avoid conflict with the names of functions
*   commonly provided by various compilers. The prefixes "C4_", "R4_" and "R8_" are used to indicate functions
*   for complex, single precision real, and double precision real arguments. 
*   For example, the sine function can be calculated by the functions C4_SIN, R4_SIN or R8_SIN.
*
*   The original, true, correct version of FN is available through NETLIB: http:  *www.netlib.org/fn/index.html.
    * Do we really need all those functions ?
    * 
    * Check what is provided by C++-0x / C99
    *
    * use overloading istead of c4 ...
    *
*
*   
*
*   Modifications by P. Hilscher (2012)
*
*   Created class FFN (Fullerton Function Library)
*   Renamed functions
*   Replaced std::complex<double> with _Complex double (should benefit vectorization)
*   Removed some common functions
*
*   @note
*   The command: from stackoverflow
*   :%s/^/foo: /
*   ...inserts foo: at the beginning of each line.
*  Reference:
*
*    Wayne Fullerton,
*    Portable Special Function Routines,
*    in Portability of Numerical Software,
*    edited by Wayne Cowell,
*    Lecture Notes in Computer Science, Volume 57,
*    Springer 1977,
*    ISBN: 978-3-540-08446-4,
*    LC: QA297.W65.
*
*    clean up with 273,280 s/\/\  * * /g ro replace comments
*
*    use separate incusion file
*
*    Set SFL Special Function Library
**/
class SFL {

  public:
  
  /**
  *   @brief evaluates the exponentially scaled Bessel function I0(X).
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument
  *    @return  the exponentially scaled Bessel function I0(X).
  *
  **/
  static double i0e( const double x);


static std::complex <float> c4_cos ( std::complex <float> z);
static std::complex <float> c4_sin ( std::complex <float> z);
static int i4_abs ( int i);
static int i4_mach ( int i);
static int i4_max ( int i1, int i2);
static int i4_min ( int i1, int i2);
static int i4_pow ( int i, int j);
static float r4_abs ( float x);
static float r4_acos ( float x);
static float r4_acosh ( float x);
static void r4_admp ( float x, float &ampl, float &phi);
static float r4_ai ( float x);
static float r4_aid ( float x);
static float r4_aide ( float x);
static float r4_aie ( float x);
static void r4_aimp ( float x, float &ampl, float &theta);
static float r4_aint ( float x);
static float r4_asin ( float x);
static float r4_asinh ( float x);
static float r4_atan ( float x);
static float r4_atan2 ( float sn, float cs);
static float r4_atanh ( float x);
static float r4_besi0 ( float x);
static float r4_besi0e ( float x);
static float r4_besi1 ( float x);
static float r4_besi1e ( float x);
static float r4_besj0 ( float x);
static float r4_besj1 ( float x);
static float r4_besk0 ( float x);
static float r4_besk0e ( float x);
static float r4_besk1 ( float x);
static float r4_besk1e ( float x);
static float *r4_beskes ( float xnu, float x, int nin);
static float *r4_besks ( float xnu, float x, int nin);
static float r4_besy0 ( float x);
static float r4_besy1 ( float x);
static float r4_beta ( float a, float b);
static float r4_betai ( float x, float pin, float qin);
static float r4_bi ( float x);
static float r4_bid ( float x);
static float r4_bide ( float x);
static float r4_bie ( float x);
static float r4_binom ( int n, int m);
static float r4_cbrt ( float x);
static float r4_chi ( float x);
static float r4_chu ( float a, float b, float x);
static float r4_chu_scaled ( float a, float b, float z);
static float r4_ci ( float x);
static float r4_cin ( float x);
static float r4_cinh ( float x);
static float r4_cos ( float x);
static float r4_cos_deg ( float x);
static float r4_cosh ( float x);
static float r4_cot ( float x);
static float r4_csevl ( float x, float a[], int n);
static float r4_dawson ( float x);
static float r4_e1 ( float x);
static float r4_ei ( float x);
static float r4_erf ( float x);
static float r4_erfc ( float x);
static float r4_exp ( float x);
static float r4_exprel ( float x);
static float r4_fac ( int n);
static float r4_gami ( float a, float x);
static float r4_gamic ( float a, float x);
static float r4_gamit ( float a, float x);
static void r4_gaml ( float &xmin, float &xmax);
static float r4_gamma ( float x);
static float r4_gamr ( float x);
static float r4_gmic ( float a, float x, float alx);
static float r4_gmit ( float a, float x, float algap1, float sgngam, float alx);
static int r4_inits ( float dos[], int nos, float eta);
static float r4_int ( float x);
static void r4_knus ( float xnu, float x, float &bknu, float &bknu1, int &iswtch);
static float r4_lbeta ( float a, float b);
static void r4_lgams ( float x, float &algam, float &sgngam);
static float r4_lgic ( float a, float x, float alx);
static float r4_lgit ( float a, float x, float algap1);
static float r4_lgmc ( float x);
static float r4_li ( float x);
static float r4_lngam ( float x);
static float r4_lnrel ( float x);
static float r4_log ( float x);
static float r4_log10 ( float x);
static float r4_mach ( int i);
static void r4_machar ( long int *ibeta, long int *it, long int *irnd, long int *ngrd,
                        long int *machep, long int *negep, long int *iexp, long int *minexp,
                        long int *maxexp, float *eps, float *epsneg, float *xmin, float *xmax);
static float r4_max ( float x, float y);
static float r4_min ( float x, float y);
static float r4_mod ( float x, float y);
static float r4_mop ( int i);
static float r4_pak ( float y, int n);
static float r4_poch ( float a, float x);
static float r4_poch1 ( float a, float x);
static float r4_power ( float a, float b);
static float r4_psi ( float x);
static float r4_rand ( float r);
static float r4_randgs ( float xmean, float sd);
static float r4_random ( float t[], int n);
static float r4_ranf ( float sw);
static float r4_ren ();
static float r4_shi ( float x);
static float r4_sign ( float x);
static float r4_si ( float x);
static void r4_sifg ( float x, float &f, float &g);
static float r4_sin ( float x);
static float r4_sin_deg ( float x);
static float r4_sinh ( float x);
static float r4_spence ( float x);
static float r4_sqrt ( float x);
static float r4_tan ( float x);
static float r4_tanh ( float x);
static void r4_upak ( float x, float &y, int &n);



   //////////////////// Double Starts /////////////////////////////////////////

  /*****************************************************************************
  *
  *  @brief i0 evaluates the Bessel function I of order 0.
  *
  *  @license  This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param  x  the argument.
  *    @return i0 the Bessel function I of order 0 of X.
  **/
  static double i0 (const double x);


  /*****************************************************************************
  *
  *    @brief evaluates the Bessel function I of order 1 of an R8 argument.
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument.
  *    @return  the Bessel function I of order 1 of X.
  **/
  static double i1 (const double x);

  /*****************************************************************************80
  *
  *    @brief evaluates the exponentially scaled Bessel function I1(X).
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument.
  *    @return  the exponentially scaled Bessel function I1(x).
  **/
  static double i1e (const double x);

  /**
  *    @brief evaluates the Bessel function J of order 0
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument.
  *    @return  the Bessel function J of order 0 of X.
  *
  **/
  static double j0 (const double x);

  /*****************************************************************************
  *
  *
  *  @brief evaluates the Bessel function J of order 1 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *
  *    @param x the argument
  *    @return  the Bessel function J of order 1 of X.
  **/
  static double j1 (const double x);

  
  /*****************************************************************************80
  *
  *   @brief evaluates the Bessel function K of order 0 of an R8 argument.
  *
  *   @license This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument.
  *    @return  the Bessel function K of order 0 of X.
  **/
  static double k0 (const double x);

  /*****************************************************************************80
  *
  *    @brief evaluates the exponentially scaled Bessel function K0(X).
  *
  *    @license This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument.
  *    @return the exponentially scaled Bessel function K0(X).
  **/
  static double k0e (const double x);

  /*****************************************************************************80
  *
  *  Purpose:
  *
  *    @brief evaluates the Bessel function K of order 1 of an R8 argument.
  *
  *    @license This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param  x the argument.
  *    @return   the Bessel function K of order 1 of X.
  **/
  static double k1 (const double x);

  /*****************************************************************************80
  *
  *    @brief evaluates the exponentially scaled Bessel function K1(X).
  *
  *    @license This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *    @param x the argument.
  *    @return  the exponentially scaled Bessel function K1(X).
  **/
  static double k1e (const double x);



static void r8_beskes ( double xnu, double x, int nin, double bke[]);
static void r8_besks ( double xnu, double x, int nin, double bk[]);


  /*****************************************************************************80
  *
  *  @brief y0 evaluates the Bessel function Y of order 0
  *
  *  @license This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *  @reference cite Wayne
  *
  *    @param x  the argument.
  *    @return   the Bessel function Y of order 0 of X
  **/
  static double y0 (const  double x);


  /*****************************************************************************80
  *
  *  @brief  y1 evaluates the Bessel function Y of order 1.
  *
  *  @license  This code is distributed under the GNU LGPL license. 
  *
  *  Author:
  *
  *    Original FORTRAN77 version by Wayne Fullerton.
  *    C++ version by John Burkardt.
  *
  *  @reference    Wayne Fullerton,
  *
  *    @param x  the argument.
  *    @return   the Bessel function Y of order 1 of X.
  **/
  static double y1 (const double x);


static double r8_beta ( double a, double b);
static double r8_betai ( double x, double pin, double qin);
static double r8_bi ( double x);
static double r8_bid ( double x);
static double r8_bide ( double x);
static double r8_bie ( double x);
static double r8_binom ( int n, int m);
static double r8_chi ( double x);
static double r8_chu ( double a, double b, double x);
static double r8_chu_scaled ( double a, double b, double z);
static double r8_ci ( double x);
static double r8_cin ( double x);
static double r8_cinh ( double x);
static double r8_cos_deg ( double x);
static double r8_cot ( double x);
static double r8_csevl ( double x, double a[], int n);
static double r8_dawson ( double x);
static double r8_e1 ( double x);
static double r8_ei ( double x);
static double r8_erf ( double x);
static double r8_erfc ( double x);
static double r8_exprel ( double x);
static double r8_fac ( int n);
static double r8_gami ( double a, double x);
static double r8_gamic ( double a, double x);
static double r8_gamit ( double a, double x);
static void r8_gaml ( double &xmin, double &xmax);
static double r8_gamma ( double x);
static double r8_gamr ( double x);
static double r8_gmic ( double a, double x, double alx);
static double r8_gmit ( double a, double x, double algap1, double sgngam, double alx);
static int r8_inits ( double dos[], int nos, double eta);
static double r8_int ( double x);
static void r8_knus ( double xnu, double x, double &bknu, double &bknu1, int &iswtch);
static double r8_lbeta ( double a, double b);
static void r8_lgams ( double x, double &algam, double &sgngam);
static double r8_lgic ( double a, double x, double alx);
static double r8_lgit ( double a, double x, double algap1);
static double r8_lgmc ( double x);
static double r8_li ( double x);
static double r8_lngam ( double x);
static double r8_lnrel ( double x);
static double r8_mach ( int i);
static void r8_machar ( long int *ibeta, long int *it, long int *irnd, long int *ngrd,
                        long int *machep, long int *negep, long int *iexp, long int *minexp,
                        long int *maxexp, double *eps, double *epsneg, double *xmin, double *xmax);
static double r8_max ( double x, double y);
static double r8_min ( double x, double y);
static double r8_mod ( double x, double y);
static double r8_mop ( int i);
static double r8_pak ( double y, int n);
static double r8_poch ( double a, double x);
static double r8_poch1 ( double a, double x);
static double r8_psi ( double x);
static double r8_ren ();
static double r8_shi ( double x);
static double r8_si ( double x);
static void r8_sifg ( double x, double &f, double &g);
static double r8_sign ( double x);
static double r8_spence ( double x);
static void r8_upak ( double x, double &y, int &n);
static void r8_admp ( double x, double &ampl, double &phi);
static double r8_ai ( double x);
static double r8_aid ( double x);
static double r8_aide ( double x);
static double r8_aie ( double x);
static void r8_aimp ( double x, double &ampl, double &theta);
static double r8_aint ( double x);
static void r8_b0mp ( double x, double &ampl, double &theta);
static void r8_b1mp ( double x, double &ampl, double &theta);



#ifdef __QUAD
  // quadrupel precision functions .....
#endif

};



#endif FFN__H_
