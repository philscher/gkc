// $Id: quadrature_gauss_1D.C 4101 2010-11-15 18:50:46Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_


#include "GaussWeights.h"

class Integrate {
Array1d points;
Array1d weights;
  
  public:

    Integrate(int len) {
      weights.resize(Range(0, len));
      points.resize(Range(0,len));

      double * points_arr = gauss_weights[len].points;
      double *weights_arr = gauss_weights[len].weights;

      for(int i = 0; i < len; i++) {
        points (i) = points_arr[i];
        weights(i) = weights_arr[i];
      }
    }

Array1d getPoints() {
    return points;

};
Array1d getWeights() {
    return weights;
};

Array3d sum(Array3d A, int n) {
    const double weight = weights(n);
    A = weight * A;
   return A; 

};

const double sum(int n) const {
    return weights(n);
};

};

#endif // __INTEGRATE_H_



