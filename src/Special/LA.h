/*
 * =====================================================================================
 *
 *       Filename: LA.cpp
 *
 *    Description: Linear Algebra Interface
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef _HELIOS_LA_H_
#define _HELIOS_LA_H_

#include "Global.h"

extern "C" {
    //int solveCycCCSymTriDiagonal(const double *b, double *x, const double *cs_alpha, const int *N);
    int solvecycccsymtridiagonal_(const double *b, double *x, const double *cs_alpha, const int *N);
}
/** 
 *  Solves a cyclic tri-diagonal system with diagonal equal to 1 and
 *  constants bands
 *  Adapted from the GNU Scientific Library (GSL) 
 * for description of method see [Engeln-Mullges + Uhlig, p. 96]
 *
 *        1         cs_alpha        0       .....  cs_alpha
 *        cs_alpha     1        cs_alpha    .....
 *        0         cs_alpha       1
 *         ....         0        cs_alpha   .....
 *        cs_alpha      ...                            1
 *
 */
  int solveCycCCSymTriDiagonal(const double *b, double *x, double cs_alpha, const int N);




#endif // __HELIOS_LA_H


