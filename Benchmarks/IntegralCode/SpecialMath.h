/* 
! Gamma/Error functions for Fortran version 2.0a
! --single-----double------quadruple---;----------defined by-----------------
!   gamma(x)   dgamma(x)   qgamma(x)   ;\int_0^\infty t^{x-1}e^{-t} dt
!   lgamma(x)  dlgamma(x)  qlgamma(x)  ;\log \Gamma(x)
!   cgamma(z)  cdgamma(z)  cqgamma(z)  ;\int_0^\infty t^{z-1}e^{-t} dt
!   clgamma(z) cdlgamma(z) cqlgamma(z) ;\log \Gamma(z)
!   erfc(x)    derfc(x)    qerfc(x)    ;2/\std::sqrt{\pi}\int_x^\infty e^{-t^2} dt
!   erf(x)     derf(x)     qerf(x)     ;2/\std::sqrt{\pi}\int_0^x e^{-t^2} dt
!   cerfc(z)   cderfc(z)   cqerfc(z)   ;2/\std::sqrt{\pi}\int_z^\infty e^{-t^2} dt
!   cerf(z)    cderf(z)    cqerf(z)    ;2/\std::sqrt{\pi}\int_0^z e^{-t^2} dt
!
 *  Author : Paul P. Hilscher (philscher@special-combo.net) 
 *  Description :
 *
 *  License : FreeBSD licesnce (see below) or GPLv3(or any later version)
 *  
 *
 *  Notes : 
 *		Copyright (c) <year>, <copyright holder>
	All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY std::expRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * License for the error (complementary) functions
    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
    You may use, copy, modify this code for any purpose and 
    without fee. You may distribute this ORIGINAL package.
 *
 *
 *
 */ 

#ifndef MATHFUNCTIONS_H__
#define MATHFUNCTIONS_H__

#include <complex>

typedef std::complex<double> cmplxd;
typedef std::complex<float>  cmplxf;


class MathFunctions {

  public:

static cmplxf erfc(cmplxf x);
static cmplxd erfc(cmplxd x);
 
 //complex quad erfc(complex quad x) {

/******************* Bessel Functions **************************/

static double i0(double x);
static double i1(double x);
  
};




#endif // MATHFUNCTIONS_H__
