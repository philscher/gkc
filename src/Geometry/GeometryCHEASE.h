/*
 * =====================================================================================
 *
 *       Filename: GeometryCHEASE.h
 *
 *    Description: Geometry using Chease numerical MHD equilibrium
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_CHEASE_H
#define GEOMETRY_CHEASE_H

#include "Geometry.h"
#include "Global.h"
#include "Setup.h"
#include "FFTSolver.h"

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


class GeometryCHEASE : public Geometry<GeometryCHEASE>
{
  double Ls;
  double shear;
 public:

   void printOn(ostream& output) const {
         output   << "Geometry  |  CHEASE (Numerical Equilibrium) Data Input : " << "mysterious file" << std::endl;
   };



  GeometryCHEASE(Setup *setup) : Geometry<GeometryCHEASE>(setup, true)  {
     

    shear    = setup->get("Geometry.Shear"   , 0.);
    chease_file = setup->get("Geometry.MHDFile", "");
    //shear = Lz / (2. * M_PI * Ls);
    //




  }

  inline const double J(const int x, const int y, const int z) { return 1.; };


  /**    Get perpendicular gradient in Fourier space. To to non-rectangular 
   *     coordiantes (shear) we have to include the non-diagonal metric component
   *     g12, g21 too
   *
  */
  inline const double k2_p(const int x_k, const int y_k, const int z) const {
      
      const double g11 = 1.;
      const double g12 = Z(z)/Ls;
      const double g22 = 1. + pow2(Z(z)/Ls);

      const double kx2 = pow2(k(Nx,Lx,x_k));
      const double ky2 = pow2(k(Ny,Ly,y_k));

      return g11 * kx2 + g22 * ky2 + 2. * g12 * sqrt(kx2 + ky2);
  };


  // for sheared magnetic fields, we have special boundary conditions, because
  // we need to take care of possible shear
  inline ShearB getYPos(const int x,const int y) {

//            // calculte shift in y, 
//            // then map by real modulo (add ghost cells)
            const double sy = - 2. * M_PI * shear * X(x);
            const int yi = ((int) (sy/dy));
            const int lDy = (NzLlD == NzGlD) ? realmod(y-3 - yi, Ny) + 3: y;
            const int uDy = (NzLuD == NzGuD) ? realmod(y-3 + yi, Ny) + 3: y;
            //const int uDy = (NzLuD == NzGuD) ? realmod(y + yi, Ly) + 3: y;
      //      std::cout << y << " " << lDy << " " << uDy << std::endl;	
//            return ShearB(y, y, 0.0  );//(yv - Y(yi-1))/dy);
            return ShearB(lDy, uDy, 0.0  );//(yv - Y(yi-1))/dy);
              
            

    check(-1, DMESG("Should not be reached here"));        
            
            
    };


private:


 int loadCheaseFile() {
    // 0) some checks
    if (mod(n_s_grid,2).eq.0) then
      call gkw_abort('N_s_grid must be odd for geom_type=chease')
    endif
    if (mod(n_s_grid,2*nperiod-1).ne.0) then
      call gkw_abort( &     'N_s_grid/(2*NPERIOD-1) has to be an integer for geom_type=chease')
      stop
    else
      npol = n_s_grid/(2*nperiod-1) !nb of points per poloidal turn for the GKW s-grid
    endif

    // used GKW file geom.f90 as reference file

    // read in number of points
   //   ! 2) READ dimensions
   //   ! number of points for psi-grid and s-grid
   //   read(igeom,*) tdum, npsi_c, tdum, ns_c
    int tdum = 
    int npsi_c =
    int tdum =
    int ns_c = 

      ! reference R and B (used for normalisation in GKW)
      ! taken to be R0EXP and B0EXO used for normalisation in CHEASE
      ! R0EXP and B0EXP close to the magnetic axis values (but not exactly)
      read(igeom,*) tdum, Rref, tdum, Bref, tdum, tdum

      ! 3) more checks
      ! check that CHEASE and GKW s-grid are compatible and that
      ! CHEASE grid is at least two times denser than GKW grid
      if (mod(ns_c,2*npol).ne.0) then
        call gkw_abort( &
   &         'NCHI in chease has to be a multiple of 2*N_s_grid/(2*NPERIOD-1)')
      else
        s_coeff=ns_c/npol ! how much the s-grid is denser in chease
      endif


      //allocate(s_indx(1:n_s_grid),stat=ierr)
      s_indx.resize(Range(1, n_s_grid));
      

      //allocate(dBdeps(1:n_s_grid),stat=ierr)
      dBdeps.resize(Range(1, n_s_grid));
      
        allocate(dBds(1:n_s_grid),stat=ierr)
      allocate(R_FS(1:n_s_grid),stat=ierr)
      allocate(Z_FS(1:n_s_grid),stat=ierr)
      allocate(dRdpsi(1:n_s_grid),stat=ierr)
      allocate(dRds(1:n_s_grid),stat=ierr)
      allocate(dZdpsi(1:n_s_grid),stat=ierr)
      allocate(dZds(1:n_s_grid),stat=ierr)
      allocate(dzetadpsi(1:n_s_grid),stat=ierr)
      allocate(dzetadchi(1:n_s_grid),stat=ierr)
      allocate(psi_dum(1:npsi_c),stat=ierr)
      allocate(psi_s_dum(1:npsi_c,1:ns_c),stat=ierr)

      ! 5) READ 1D arrays

      ! find the index for the radial grid
      switch (eps_type) {
          // use eps to select FS
        case(1) :
          // psi-grid, s-grid, Rgeom (not used)
          read(igeom,*) tdum, (dum,i=1,npsi_c) !psi
          read(igeom,*) tdum, (dum,i=1,ns_c) !s
          read(igeom,*) tdum, (dum,i=1,npsi_c) !Rgeom
          // amin=(Rmax-Rmin)/2
          read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) !amin
          // find the index for the psi grid
          psi_indx=1
          do while (abs(psi_dum(psi_indx)/Rref-eps)>0.008*psi_dum(npsi_c)/Rref .and. psi_indx.LT.npsi_c)
            psi_indx = psi_indx+1
          end do

        case (2)! use rho_psi to select FS
          // psi-grid
          read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) !psi
          // find the index for the psi grid
          psi_indx=1
          do while (abs(sqrt(psi_dum(psi_indx)/psi_dum(npsi_c))-eps)>0.008 .and. psi_indx.LT.npsi_c)
            psi_indx = psi_indx+1
          end do

          // s-grid, Rgeom (not used)
          read(igeom,*) tdum, (dum,i=1,ns_c) !s
          read(igeom,*) tdum, (dum,i=1,npsi_c) !Rgeom

          // amin=(Rmax-Rmin)/2
          read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) !amin

          default : check(-1, DMESG("Value not allowed for eps_type"));
      }
      if (psi_indx .GE. npsi_c) check(-1, DMESG("Radial grid too coarse in CHEASE: increase NPSI"));

      ! value of eps=amin/R (used in mode.F90 for kxspace)
      eps = psi_dum(psi_indx)/Rref

      ! depsdpsi
      ! -> to go from chease (psi,zeta,s) to GKW (eps,zeta,s)
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! damindpsi
      depsdpsi = psi_dum(psi_indx)/Rref

      ! Bmax, Bmin, q, shat
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! Bmax
      bmax = psi_dum(psi_indx)/Bref
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! Bmin
      bmin = psi_dum(psi_indx)/Bref
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! q
      q = psi_dum(psi_indx)
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! dqdpsi
      dqdpsi = psi_dum(psi_indx)
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! dqdpsi_chk
      dqdpsi_dum = psi_dum(psi_indx)
      shat = eps * dqdpsi_dum / q / depsdpsi

      ! p and dpdeps (not normalised)
      ! -> to be used for beta and beta_prime
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! p
      p = psi_dum(psi_indx)
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! dpdpsi
      dpdeps = psi_dum(psi_indx)/depsdpsi

      ! jacobian J_psi_zeta_s
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c) ! jacobian
      jac=psi_dum(psi_indx)

      ! not used (don't change the order!)
      read(igeom,*) tdum, (dum,i=1,npsi_c) ! djacdpsi

      !chease F (not the same as GKW F tensor which is ffun)
      read(igeom,*) tdum, (psi_dum(i),i=1,npsi_c)
      F=psi_dum(psi_indx)

      !not used
      read(igeom,*) tdum, (dum,i=1,npsi_c) ! dFdpsi

      ! 6) READ 2D arrays
      ! build the index correspondance for the s-grid
      ! CHEASE s-array goes from 0 to 1 (LFS midplane, counterclockwise)
      do i = 1, npol
        do j = 1, 2*nperiod - 1
          dum=s_coeff*(i-1)+s_coeff/2+1
          if (dum.gt.ns_c/2) then
            s_indx(i+(j-1)*npol)=dum-ns_c/2
          else
            s_indx(i+(j-1)*npol)=dum+ns_c/2
          endif
        end do
      end do

      ! metric_G elements, (psi,zeta,s) coordinates
      ! g_psi_psi
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !g11
      do i = 1, n_s_grid
        metric_G(i,1,1)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! g_psi_zeta
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !g12
      do i = 1, n_s_grid
        metric_G(i,1,2)=psi_s_dum(psi_indx,s_indx(i))
        metric_G(i,2,1)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! g_psi_s
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !g13
      do i = 1, n_s_grid
        metric_G(i,1,3)=psi_s_dum(psi_indx,s_indx(i))
        metric_G(i,3,1)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! g_zeta_zeta
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !g22
      do i = 1, n_s_grid
        metric_G(i,2,2)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! g_zeta_s
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !g23
      do i = 1, n_s_grid
        metric_G(i,2,3)=psi_s_dum(psi_indx,s_indx(i))
        metric_G(i,3,2)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! g_s_s
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !g33
      do i = 1, n_s_grid
        metric_G(i,3,3)=psi_s_dum(psi_indx,s_indx(i))
      end do


      ! magnetic field
      ! norm(B)
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) !B
      do i = 1, n_s_grid
        bn_G(i)=psi_s_dum(psi_indx,s_indx(i))/Bref
      end do
      ! dBdeps (not normalized)
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dBdpsi
      do i = 1, n_s_grid
        dBdeps(i)=psi_s_dum(psi_indx,s_indx(i))/depsdpsi
      end do
      ! dBds (not normalized)
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dBds
      do i = 1, n_s_grid
        dBds(i)=psi_s_dum(psi_indx,s_indx(i))
      end do

      ! R and Z
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! R
      do i = 1, n_s_grid
        R_FS(i)=psi_s_dum(psi_indx,s_indx(i))
      end do
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! Z
      do i = 1, n_s_grid
        Z_FS(i)=psi_s_dum(psi_indx,s_indx(i))
      end do

      ! R and Z gradients
      ! dRdpsi
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dRdpsi
      do i = 1, n_s_grid
        dRdpsi(i)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! dRds
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dRds
      do i = 1, n_s_grid
        dRds(i)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! dZdpsi
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dZdpsi
      do i = 1, n_s_grid
        dZdpsi(i)=psi_s_dum(psi_indx,s_indx(i))
      end do
      ! dZds
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dZds
      do i = 1, n_s_grid
        dZds(i)=psi_s_dum(psi_indx,s_indx(i))
      end do

      ! poloidal angle
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! theta
      do i = 1, n_s_grid
        dum=psi_s_dum(psi_indx,s_indx(i))
        if (dum.gt.pi) then
          dum = dum - 2E0*pi
        endif
        pol_angle(i) = dum + 2E0*pi*(floor(real(i-1)/real(npol))-real(nperiod)+1)
      end do

      ! not used
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((dum,i=1, npsi_c),j=1, ns_c) ! Ah
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((dum,i=1, npsi_c),j=1, ns_c) ! dAhdpsi

      ! dzetadpsi
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dzetadpsi
      do i = 1, n_s_grid
        dzetadpsi(i)=psi_s_dum(psi_indx,s_indx(i))
      end do

      ! dzetadchi
      read(igeom,'(A)') tdum
      read(igeom,'(1P5E20.10)') ((psi_s_dum(i,j),i=1, npsi_c),j=1, ns_c) ! dzetadchi
      do i = 1, n_s_grid
        dzetadchi(i)=psi_s_dum(psi_indx,s_indx(i))
      end do

      close(igeom)


      ! 7) Compute the missing elements

      ! Correction to the metric_G elements involving dzetadpsi
      ! because dzetadpsi is not periodic in s:
      ! dzetadpsi(1) = dqdpsi_dum
      ! dzetadpsi(s) = dzetadpsi(s_0) + dum * dqdpsi
      !    with  s=s_0+ dum  and  0<=s_0<1
      do i = 1, n_s_grid
        dum = int(i/npol) - nperiod
        if (mod(i,npol).GT.real(npol/2)) then
          dum = dum + 1
        endif
        metric_G(i,1,2) = metric_G(i,1,2) + &
                      & dum * dqdpsi_dum * metric_G(i,1,1)
        metric_G(i,2,1) = metric_G(i,2,1) + &
                      & dum * dqdpsi_dum * metric_G(i,1,1)
        metric_G(i,2,3) = metric_G(i,2,3) + &
                      & dum * dqdpsi_dum * metric_G(i,1,3)
        metric_G(i,3,2) = metric_G(i,3,2) + &
                      & dum * dqdpsi_dum * metric_G(i,1,3)
        metric_G(i,2,2) = metric_G(i,2,2) + (dum * dqdpsi_dum)**2 * metric_G(i,1,1) + &
                      & 2 * dum * dqdpsi_dum * dzetadpsi(i) * metric_G(i,1,1) + &
                      & 4 * pi * dum * dqdpsi_dum * dzetadchi(i) * metric_G(i,1,3)
      end do

      ! beta and beta prime (does not have the bn_G dependence)
      beta_real = 2 * 1.25663706144E-6 * p / Bref**2
      beta_prime_real = 2 * 1.25663706144E-6 * dpdeps / Bref**2

      ! calculate the array for the parallel derivative (ffun)
      do i = 1, n_s_grid
        ffun(i) = 2*pi*signJ*Rref/bn_G(i)/Bref/jac
      end do

      ! the function connected with the trapping terms (gfun)
      do i = 1, n_s_grid
        gfun(i) = ffun(i)*dBds(i)/bn_G(i)/Bref
      end do

      ! calculate the array connected with the ExB velocity (efun)
      do i = 1, n_s_grid

        ! the diagonal components are zero
        efun(i,1,1) = 0.
        efun(i,2,2) = 0.
        efun(i,3,3) = 0.

        ! the psi zeta component
        efun(i,1,2) =  signJ*pi*Rref**2/bn_G(i)**2/Bref*  &
                    & (metric_G(i,1,1)*metric_G(i,2,2) -  &
                    &  metric_G(i,1,2)**2)*depsdpsi

        ! the psi s component
        efun(i,1,3) =  signJ*pi*Rref**2/bn_G(i)**2/Bref*  &
                    & (metric_G(i,1,1)*metric_G(i,2,3) -  &
                    &  metric_G(i,1,2)*metric_G(i,1,3))*depsdpsi*signB*signJ

        ! the zeta s component
        efun(i,2,3) =  signJ* pi*Rref**2/bn_G(i)**2/Bref* &
                    & (metric_G(i,1,2)*metric_G(i,2,3) -  &
                    &  metric_G(i,2,2)*metric_G(i,1,3))

        ! the other components are anti-symmetric_G
        efun(i,2,1) = - efun(i,1,2)
        efun(i,3,1) = - efun(i,1,3)
        efun(i,3,2) = - efun(i,2,3)

      end do

      ! calculate the curvature function  (dfun)
      do i = 1, n_s_grid

        ! the psi component
        dfun(i,1) = -2/bn_G(i)/Bref*efun(i,1,3)*dBds(i)

        ! the zeta component
        dfun(i,2) = -2/bn_G(i)/Bref*(efun(i,2,1)*dBdeps(i) + &
                  &  efun(i,2,3)*dBds(i))

        ! the s component
        dfun(i,3) = -2/bn_G(i)/Bref*efun(i,3,1)*dBdeps(i)

      end do

      ! calculate the array connected with the coriolis drift (hfun)
      ! Omega/vthref=Cte assumed (rigid body)
      ! As VCOR=vtor(R=Rref)/vthref, Omega/vthref=VCOR/Rref
      ! term VCOR=vtor(R=Rref)/vthref not included here (added in linear_terms)
      do i = 1, n_s_grid

        ! the psi component
        hfun(i,1) = - signB*Rref**2/bn_G(i)*(dZdpsi(i)*metric_G(i,1,1) + &
                  & dZds(i)*metric_G(i,3,1))*depsdpsi / Rref

        ! the zeta component
        hfun(i,2) = - signB*Rref**2/bn_G(i)* &
                  & (dZdpsi(i)*metric_G(i,1,2)*signB*signJ + &
                  & dZds(i)*metric_G(i,3,2)*signB*signJ) / Rref

        ! the s component
        hfun(i,3) = - signB*Rref**2/bn_G(i)*(dZdpsi(i)*metric_G(i,1,3) + &
                  & dZds(i)*metric_G(i,3,3) - dZds(i)*ffun(i)**2/Rref**2) &
                  & / Rref

      end do

      ! calculate the array connected with centrifugal drift (ifun)
      ! does not include the (Rref*Omega/vthref)**2 term
      do i = 1, n_s_grid

        ! the psi component
        ifun(i,1) = 2*R_FS(i)*(efun(i,1,3)*dRds(i))/ Rref**2

        ! the zeta component
        ifun(i,2) = 2*R_FS(i)*(efun(i,2,1)*dRdpsi(i)/depsdpsi+ &
                  & efun(i,2,3)*dRds(i))/ Rref**2
        ! the s component
        ifun(i,3) = 2*R_FS(i)*(efun(i,3,1)*dRdpsi(i)/depsdpsi)/ Rref**2

      end do

      ! Calculate the normfactor for k_zeta
      kthnorm = Rref * sqrt(metric_G((n_s_grid+1)/2,2,2))

      ! metric_G elements for the (eps,zeta,s) coordinates
      ! all multiplied by Rref**2 for normalisation
      do i=  1, n_s_grid
        metric_G(i,1,1) = metric_G(i,1,1) * depsdpsi**2 * Rref**2
        metric_G(i,1,2) = signB*signJ*metric_G(i,1,2) * depsdpsi * Rref**2
        metric_G(i,2,1) = signB*signJ*metric_G(i,2,1) * depsdpsi * Rref**2
        metric_G(i,1,3) = metric_G(i,1,3) * depsdpsi * Rref**2
        metric_G(i,3,1) = metric_G(i,3,1) * depsdpsi * Rref**2
        metric_G(i,2,2) = metric_G(i,2,2) * Rref**2
        metric_G(i,2,3) = signB*signJ*metric_G(i,2,3) * Rref**2
        metric_G(i,3,2) = signB*signJ*metric_G(i,3,2) * Rref**2
        metric_G(i,3,3) = metric_G(i,3,3) * Rref**2
      end do

      !Compute the toroidal and poloidal fractions of the field
      !CHECK IF THIS IS CORRECT, did not have details of hamada.dat normalisations
      !Done by trial / error, seems right..
      do i=  1, n_s_grid
        !bt_frac should always be positive
        bt_frac(i)=signB*F/(R_FS(i)*Bref*bn_G(i))
        !bt_frac should be negative for negative signJ
        bp_frac(i)=signJ*sqrt(1-bt_frac(i)**2)

        !Some checks
        if(bt_frac(i)<0..or.bt_frac(i)>1) then
          !call gkw_abort('error in bt_frac')
        end if

        if(signJ*bp_frac(i)<0.or.bp_frac(i)>1) then
           !call gkw_abort('error in bt_frac')
        end if
      end do

      !NOT YET checked, for now revert to s-alpha usage!!!!!
      bt_frac(:)=1.
      bp_frac(:)=signJ*eps/q

    !end if !root processor

    call mpibcast_real(Bref,            1)
    call mpibcast_real(Rref,            1)
    call mpibcast_real(q,               1)
    call mpibcast_real(shat,            1)
    call mpibcast_real(eps,             1)
    call mpibcast_real(bmin,            1)
    call mpibcast_real(bmax,            1)
    call mpibcast_real(kthnorm,         1)
    call mpibcast_real(beta_real,       1)
    call mpibcast_real(beta_prime_real, 1)
    call mpibcast_real(bn_G,            n_s_grid)
    call mpibcast_real(pol_angle,       n_s_grid)
    call mpibcast_real(metric_G,      9*n_s_grid)
    call mpibcast_real(dfun,          3*n_s_grid)
    call mpibcast_real(efun,          9*n_s_grid)
    call mpibcast_real(ffun,            n_s_grid)
    call mpibcast_real(gfun,            n_s_grid)
    call mpibcast_real(hfun,          3*n_s_grid)
    call mpibcast_real(ifun,          3*n_s_grid)
    call mpibcast_real(bt_frac,          n_s_grid)
    call mpibcast_real(bp_frac,          n_s_grid)

  case default
    ! No known option specified
    write(*,*)'You specified ',geom_type, &
    & 'for geom_type in the namelist GEOM'
    write(*,*)'Only known options are: s-alpha, chease'
    stop

    }
#endif // GEOMETRY_H
