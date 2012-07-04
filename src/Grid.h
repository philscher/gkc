/*
 * =====================================================================================
 *
 *       Filename: Grid.h
 *
 *    Description: Grid definitions incuding boundary
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GRID_H_
#define __GRID_H_


#include "config.h"
#include "Setup.h"
#include "Parallel.h"
#include "Integration.h"
#include "Global.h"

#include "FileIO.h"
//! Grid - Class setting basic properties and size of the computational domain
/*! The compuational domain consist basically of..
*   \image html grid.png
    For the Parallized these boudaries are transfered between the
*   CPUs using MPI. We use 2-D decomposition in X and Y. Because the
*
*
*/

class Grid : public IfaceHelios {

public:
Array4d  SendYu, SendXu, SendYl, SendXl; 
Array4d  RecvYu, RecvXu, RecvYl, RecvXl;

Integration *Ipol_M;

// GhostCell
int NxGC, NyGC, NzGC, NvGC;

Range RxGB, RyGB, RzGB, RvGB, RmGB, RsGB; 
Range RxGD, RyGD, RzGD, RvGD, RmGD, RsGD; 


// Local Grid size (for Domain & Boundary)

      int NxGD, NyGD, NkyGD, NzGD, NvGD, NmGD, NsGD;


  Grid(Setup *setup, Parallel *parallel, FileIO *fileIO);
 ~Grid();
protected:
     virtual void printOn(ostream &output) const;


public: 
    void initDataOutput(FileIO *fileIO);
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

     int getGlobalSize() const { return Nx * Nky * Nz * Nv * Nm * Ns; };
     int getLocalSize() const { return NxLD * NkyLD * NzLD * NvLD * NmLD * NsLD;};
    
};



#endif // __GRID_H_

