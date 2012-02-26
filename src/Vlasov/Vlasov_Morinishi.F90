!=======================================================================
!     subroutine morinishi
!
!          (1)calculate the time derivatve by 4th order Morinishi scheme
!          (2)time-integrate the distribution function by 4th order
!             Runge-Kutta scheme
!=======================================================================


#include "config.h"

module Morinishi
   use ISO_C_BINDING
 implicit none
      integer(C_INT), bind(C, name="NxLlD") :: NxLlD
      integer(C_INT), bind(C, name="NxLuD") :: NxLuD
      integer(C_INT), bind(C, name="NyLlD") :: NyLlD
      integer(C_INT), bind(C, name="NyLuD") :: NyLuD
      integer(C_INT), bind(C, name="NzLlD") :: NzLlD
      integer(C_INT), bind(C, name="NzLuD") :: NzLuD
      integer(C_INT), bind(C, name="NvLlD") :: NvLlD
      integer(C_INT), bind(C, name="NvLuD") :: NvLuD
      integer(C_INT), bind(C, name="NmLlD") :: NmLlD
      integer(C_INT), bind(C, name="NmLuD") :: NmLuD
      integer(C_INT), bind(C, name="NsLlD") :: NsLlD
      integer(C_INT), bind(C, name="NsLuD") :: NsLuD
      
      integer(C_INT), bind(C, name="NxLlB") :: NxLlB
      integer(C_INT), bind(C, name="NxLuB") :: NxLuB
      integer(C_INT), bind(C, name="NyLlB") :: NyLlB
      integer(C_INT), bind(C, name="NyLuB") :: NyLuB
      integer(C_INT), bind(C, name="NzLlB") :: NzLlB
      integer(C_INT), bind(C, name="NzLuB") :: NzLuB
      integer(C_INT), bind(C, name="NvLlB") :: NvLlB
      integer(C_INT), bind(C, name="NvLuB") :: NvLuB
      integer(C_INT), bind(C, name="NmLlB") :: NmLlB
      integer(C_INT), bind(C, name="NmLuB") :: NmLuB
      integer(C_INT), bind(C, name="NsLlB") :: NsLlB
      integer(C_INT), bind(C, name="NsLuB") :: NsLuB

      integer(C_INT), bind(C, name="NzGlB") :: NzGlB
      integer(C_INT), bind(C, name="NzGuB") :: NzGuB
      integer(C_INT), bind(C, name="NvGlB") :: NvGlB
      integer(C_INT), bind(C, name="NvGuB") :: NvGuB
      integer(C_INT), bind(C, name="NmGlB") :: NmGlB
      integer(C_INT), bind(C, name="NmGuB") :: NmGuB
      
      logical(C_BOOL), bind(C, name="do_gyro") :: do_gyro

      real (C_DOUBLE),  bind(C) :: dx, dy, dz, dv, dm
      
      double precision :: c_speed = 2.99792458e10

      type, BIND(C) :: species2
            real (C_DOUBLE) :: charge
            real (C_DOUBLE) :: mass
            real (C_DOUBLE) :: collision
            real (C_DOUBLE) :: w_n
            real (C_DOUBLE) :: w_T
            logical(C_BOOL) :: gyroAverage
            real (C_DOUBLE) :: T0
            real (C_DOUBLE) :: n0
            real (C_DOUBLE) :: scale_v
            real (C_DOUBLE) :: scale_n
            real (C_DOUBLE) :: sigma
            real (C_DOUBLE) :: alpha
            character(64)   :: name  
            character(64)   :: T_name  
            character(64)   :: n_name  
            ! Note : Its a pointer
            character(40) :: T
            character(40) :: n
    !std::string name;
      END TYPE species2

!      INTERFACE 
!        SUBROUTINE solveYTriDiagonal(b,x,cs_alpha, N) BIND(C,NAME="solveCycCCSymTriDiagonal")
!            USE, INTRINSIC :: iso_c_binding, ONLY : c_double, c_int
!            REAL (c_double)  :: b(:),x(:),cs_alpha
!            INTEGER(c_int)  :: N
!        END SUBROUTINE solveYTriDiagonal
!      END INTERFACE


contains

subroutine Fortran_Vlasov_EM_NonLinear(f0, f, fs, fss, ft, H, Xi, phi &
, Ap, dt, rk_step, sp, Vel, Mu, B0, beta, eps, max_Xi) bind(C, name="Fortran_Vlasov_EM_NonLinear")

      implicit none

      integer :: rk_step
      double precision  ::  q, mass, sub = 1.d0/2.d0, B0, w_n, w_T, scale_v, scale_n, sigma, alpha, &
                            beta, eps, dfs_dv, dH_dz
      double precision  ::  Vel(NvGlB:NvGuB),  Mu(NmGlB:NmGuB), T, N, &
      phi(NxLlB-2:NxLuB+2, NyLlB-2:NyLuB+2, NzLlB:NzLuB, NmLlD:NmLuD, NsLlD:NsLuD),                       &
      Ap (NxLlB-2:NxLuB+2, NyLlB-2:NyLuB+2, NzLlB:NzLuB, NmLlD:NmLuD, NsLlD:NsLuD),                       &
      Xi (NxLlB-2:NxLuB+2, NyLlB-2:NyLuB+2, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
      H(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
      f(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
      f0(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),          &
      fs(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),          &
      fss(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),         &
      ft(NxLlD:NxLuD, NyLlD:NyLuD, NzLlD:NzLuD, NvLlD:NvLuD, NmLlD:NmLuD, NsLlD:NsLuD), fts, dt

      double precision :: Xi_y, Xi_y_xm1, Xi_y_xm2, Xi_y_xp1, Xi_y_xp2, Xi_x,&
      Xi_x_ym1, Xi_x_ym2, Xi_x_yp1, Xi_x_yp2, Xi_z
      integer :: v, z, y, x, s, m, si =1
    
      double precision :: max_Xi(0:2)

      type(species2) :: sp(0:8)
      
      max_Xi(:) = 0.
   
      if(do_gyro .eqv. .TRUE.) sub = 3.0d0/2.0d0
     
     do s = NsLlD, NsLuD
       ! charge over mass
       q     = sp(s)%charge
       mass  = sp(s)%mass
       w_n   = sp(s)%w_n
       w_T   = sp(s)%w_T
       T     = sp(s)%T0
       n     = sp(s)%n0
       scale_v    = sp(s)%scale_v
       scale_n    = sp(s)%scale_n
       sigma    = sp(s)%sigma
       alpha    = sp(s)%alpha
    
      
        do m = NmLlB, NmLuB
            do v = NvLlB, NvLuB
                do  z = NzLlB, NzLuB
                    do  y = NyLlB-2, NyLuB+2
                        do  x = NxLlB-2, NxLuB+2
                            Xi(x,y,z,v,m,s) = phi(x,y,z,m,s) !- alpha * eps * beta * Vel(v)*Ap(x,y,z,m,s)
                        enddo
                    enddo
                    
                    do  y = NyLlB, NyLuB
                        do  x = NxLlB, NxLuB
                            H(x,y,z,v,m,s)  = fs(x,y,z,v,m,s) + sigma * Xi(x,y,z,v,m,s) * f0(x,y,z,v,m,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
       


     do s = NsLlD, NsLuD
       ! charge over mass
       q     = sp(s)%charge
       mass  = sp(s)%mass
       w_n   = sp(s)%w_n
       w_T   = sp(s)%w_T
       T     = sp(s)%T0
       n     = sp(s)%n0
       scale_v    = sp(s)%scale_v
       scale_n    = sp(s)%scale_n

        do m = NmLlD, NmLuD
            do v = NvLlD, NvLuD
                do  z = NzLlD, NzLuD
                    do  y = NyLlD, NyLuD
                        do  x = NxLlD, NxLuD


                      ! calculate here the derivatives otherwise we would have to allocate whole array, how (in)efficient is this ?
                      ! calculate \frac{\Xe}{x} 
                      ! probably we can improve this by using cycling shift ? 
                      Xi_y     = (2.d0/ 3.d0*(Xi(x  , y+1, z, v, m, s) - Xi(x  , y-1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x  , y+2, z, v, m, s) - Xi(x  , y-2, z, v, m, s)))/dy    
                      Xi_y_xm1 = (2.d0/ 3.d0*(Xi(x-1, y+1, z, v, m, s) - Xi(x-1, y-1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x-1, y+2, z, v, m, s) - Xi(x-1, y-2, z, v, m, s)))/dy    
                      Xi_y_xm2 = (2.d0/ 3.d0*(Xi(x-2, y+1, z, v, m, s) - Xi(x-2, y-1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x-2, y+2, z, v, m, s) - Xi(x-2, y-2, z, v, m, s)))/dy    
                      Xi_y_xp1 = (2.d0/ 3.d0*(Xi(x+1, y+1, z, v, m, s) - Xi(x+1, y-1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x+1, y+2, z, v, m, s) - Xi(x+1, y-2, z, v, m, s)))/dy    
                      Xi_y_xp2 = (2.d0/ 3.d0*(Xi(x+2, y+1, z, v, m, s) - Xi(x+2, y-1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x+2, y+2, z, v, m, s) - Xi(x+2, y-2, z, v, m, s)))/dy    

                      Xi_x     = (2.d0/ 3.d0*(Xi(x+1, y  , z, v, m, s) - Xi(x-1, y  , z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x+2, y  , z, v, m, s) - Xi(x-2, y  , z, v, m, s)))/dx     
                      Xi_x_ym1 = (2.d0/ 3.d0*(Xi(x+1, y-1, z, v, m, s) - Xi(x-1, y-1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x+2, y-1, z, v, m, s) - Xi(x-2, y-1, z, v, m, s)))/dx     
                      Xi_x_ym2 = (2.d0/ 3.d0*(Xi(x+1, y-2, z, v, m, s) - Xi(x-1, y-2, z, v, m, s))   &
                                 -1.d0/12.d0*(Xi(x+2, y-2, z, v, m, s) - Xi(x-2, y-2, z, v, m, s)))/dx     
                      Xi_x_yp1 = (2.d0/ 3.d0*(Xi(x+1, y+1, z, v, m, s) - Xi(x-1, y+1, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x+2, y+1, z, v, m, s) - Xi(x-2, y+1, z, v, m, s)))/dx     
                      Xi_x_yp2 = (2.d0/ 3.d0*(Xi(x+1, y+2, z, v, m, s) - Xi(x-1, y+2, z, v, m, s))  &
                                 -1.d0/12.d0*(Xi(x+2, y+2, z, v, m, s) - Xi(x-2, y+2, z, v, m, s)))/dx     

                      Xi_z     = (2.d0/ 3.d0*(Xi(x, y, z+1, v, m, s) - Xi(x, y, z-1, v, m, s)) &
                                 -1.d0/12.d0*(Xi(x, y, z+2, v, m, s) - Xi(x, y, z-2, v, m, s)))/dz 
   
                      dfs_dv   = (2.d0/ 3.d0*(fs(x, y, z, v+1, m, s) - fs(x, y, z, v-1, m, s))               &
                                 -1.d0/12.d0*(fs(x, y, z, v+2, m, s) - fs(x, y, z, v-2, m, s)))/dv

                       dH_dz   = (4.d0/ 3.d0*(H (x, y, z+1, v, m, s) - H (x, y, z-1, v, m, s))               &
                                 -1.d0/12.d0*(H (x, y, z+2, v, m, s) - H (x, y, z-2, v, m, s)))/dz


                        if(rk_step == 4) then 
                            max_Xi(0) = max(max_Xi(0), abs(Xi_x))
                            max_Xi(1) = max(max_Xi(1), abs(Xi_y))
                            max_Xi(2) = max(max_Xi(2), abs(Xi_z))
                        endif
        ! Finally the Vlasov equation
        fts = &
      ! calculate the time derivatve      
           !    Modfied non-linear terms 
            ( 1.d0/ 3.d0*((Xi_y+Xi_y_xp1)*fs(x+1, y, z, v, m, s) - (Xi_y+Xi_y_xm1)*fs(x-1, y, z, v, m, s))      &
            - 1.d0/24.d0*((Xi_y+Xi_y_xp2)*fs(x+2, y, z, v, m, s) - (Xi_y+Xi_y_xm2)*fs(x-2, y, z, v, m, s)))/dx     &
          - ( 1.d0/3.d0*((Xi_x+Xi_x_yp1)*fs(x, y+1, z, v, m, s) - (Xi_x+Xi_x_ym1)*fs(x, y-1, z, v, m, s))      &
            - 1.d0/24.d0*((Xi_x+Xi_x_yp2)*fs(x, y+2, z, v, m, s) - (Xi_x+Xi_x_ym2)*fs(x, y-2, z, v, m, s)))/dy     &
          ! Term 3
          - Vel(v)*scale_v*dH_dz   &
          ! Term 4
          - Xi_y * (w_n + w_T * (((Vel(v)**2.e0+Mu(m)* B0)  - sub)))  * f0(x, y, z, v, m, s ) &
          + sigma/mass/scale_v**2 * Xi_z* dfs_dv

        !   time-integrate the distribution function     
        if(rk_step == 1) then
            ft(x, y, z, v, m, s) = fts
        else if((rk_step == 2) .or. (rk_step == 3)) then
            ft(x, y, z, v, m, s) = ft(x, y, z, v, m, s) + 2.d0*fts
        else    
            fts = ft(x, y, z, v, m, s) + fts
        end if
        
        fss(x, y, z, v, m, s) = f(x, y, z, v, m, s) + fts*dt
          
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    return
end subroutine Fortran_Vlasov_EM_NonLinear


subroutine morinishi_sa(Vel, Zv, f0, f, fs, fss, ft, Ex, Ey, Ez, dt, rk_step, sp, Mu, B0, shear, R0, q0, &
epsilon_t, use_terms) bind(C, name="Morinishi_SA")

      implicit none

      integer :: rk_step, use_terms
      double precision  ::  q, mass, sub = 1.d0/2.d0, B0, w_n, w_T, scale_v, scale_n           
      double precision :: R0, q0, epsilon_t, shear
      double precision  ::  Vel(NvGlB:NvGuB),  Mu(NmGlB:NmGuB), Zv(NzGlB:NzGuB), T, N, &
      Ex(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NmLlB:NmLuB, NsLlB:NsLuB),                       &
      Ey(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NmLlB:NmLuB, NsLlB:NsLuB),                       &
      Ez(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NmLlB:NmLuB, NsLlB:NsLuB),                       &
      f(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
      f0(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),          &
      fs(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),          &
      fss(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),         &
      ft(NxLlD:NxLuD, NyLlD:NyLuD, NzLlD:NzLuD, NvLlD:NvLuD, NmLlD:NmLuD, NsLlD:NsLuD), fts, dt

      double precision :: Kx, Ky, Lo = 1., sigma=1., alpha_MHD = 0., cs=1., alpha_s=1.
      integer :: v, z, y, x, s, m, si =1

      logical :: geo_term = .FALSE., nonlin_term = .TRUE., lin_term = .TRUE., ptrap_term = .FALSE., vpar_term = .TRUE.
     
      ! set terms
      
      
      
      type(species2) :: sp(0:8)
      

!      lin_term    = IAND(use_terms,  1 )
!      nonlin_term = IAND(use_terms,  2 )
!      geo_term    = IAND(use_terms,  4 )
!      ptrap_term  = IAND(use_terms,  8 )
!      vpar_term   = IAND(use_terms,  16 )

      if(do_gyro .eqv. .TRUE.) sub = 3.0d0/2.0d0

      do s = NsLlD, NsLuD
       ! charge over mass
       q     = sp(s)%charge
       mass  = sp(s)%mass
       w_n   = sp(s)%w_n
       w_T   = sp(s)%w_T
       T     = sp(s)%T0
       n     = sp(s)%n0
       scale_v    = sp(s)%scale_v
       scale_n    = sp(s)%scale_n
       sigma      = q/T

       alpha_s =scale_v !Lo / (q0/R0) * scale_v / cs
       if(do_gyro .eqv. .TRUE.) si = s
        do m = NmLlD, NmLuD
            do v = NvLlD, NvLuD
                do  z = NzLlD, NzLuD
          !          Kx = - 2. * Lo/R0 * sin(Zv(z))
                    ! includes Shafranov shift, see Danners PhD Thesis, Eq. (2.37)
           !         Ky = - 2.* Lo/R0 * (cos(Zv(z)) + (shear * Zv(z) - alpha_MHD * sin(Zv(z)))* sin(Zv(z)))

! shear
                       Kx = 0.
                       Ky = shear * Zv(z)

                    do  y = NyLlD, NyLuD
                        do  x = NxLlD, NxLuD

      ! calculate the time derivatve      
      fts = 0.
      ! Non-Linear ExB Terms + vp-nonlinearity
      if(lin_term) then
            fts = fts  &
          ! Driving from Local assumption
          - Ey(x, y, z, m,s) *(w_n + w_T * (((Vel(v)**2.e0+Mu(m)* B0)  - sub)))  * f0(x, y, z, v, m, s ) &
         ! Velocity 
          - alpha_s * Vel(v) * &
            (   4.d0/3.d0 * (fs(x, y, z+1, v, m, s) - fs(x, y, z-1, v, m, s))/(2.d0*dz)               &
              - 1.d0/3.d0 * (fs(x, y, z+2, v, m, s) - fs(x, y, z-2, v, m, s))/(4.d0*dz)               &
            )               &
          ! Driving from local assumption
          - alpha_s*Vel(v) * sigma  * Ez(x, y, z, m,s)* f0(x, y, z, v, m, s)  
      endif


      if(nonlin_term) then
            fts = fts + &
         ! Term 2
            4.d0/3.d0                                      &
          *(Ey(x+1, y, z, m,s)* fs(x+1, y, z, v, m, s) - Ey(x-1, y, z, m,s)*fs(x-1, y, z, v, m, s)            &
          + Ey(x  , y, z, m,s)*(fs(x+1, y, z, v, m, s) - fs(x-1, y, z, v, m, s)))/(4.d0*dx)              &
          - 1.d0/3.d0                                       &
          *(Ey(x+2, y, z, m,s)* fs(x+2, y, z, v, m, s) - Ey(x-2, y, z, m,s)*fs(x-2, y, z, v, m, s)            &
          + Ey(x  , y, z, m,s)*(fs(x+2, y, z, v, m, s) - fs(x-2, y, z, v, m, s)))/(8.d0*dx)              &
          ! Term 1
          - (4.d0/3.d0                                       &
          *(Ex(x, y+1, z, m,s)*fs(x, y+1, z, v, m, s) - Ex(x, y-1, z, m,s)*fs(x, y-1, z, v, m, s)            &
         +  Ex(x, y  , z, m,s)*(fs(x, y+1, z, v, m, s)  - fs(x, y-1, z, v, m, s)))/(4.d0*dy)              &
          - 1.d0/3.d0                                       &
          *(Ex(x, y+2, z, m,s)* fs(x, y+2, z, v, m, s) - Ex(x, y-2, z, m,s)*fs(x, y-2, z, v, m, s)            &
         +  Ex(x, y  , z, m,s)*(fs(x, y+2, z, v, m, s) - fs(x, y-2, z, v, m, s)))/(8.d0*dy))              
       endif
       
       ! Geometric terms
       if (geo_term) then
         !write(*,*) " geom"
            fts = fts  - (Mu(m)*B0+ 2.*Vel(v)**2)/(2.*sigma*B0) * ( &
         ! Term 2
           Kx *  & 
          ( 4.d0/3.d0 * (fs(x+1, y, z, v, m, s) - fs(x-1, y, z, v, m, s))/(2.d0*dx)       &
          - 1.d0/3.d0 * (fs(x+2, y, z, v, m, s) - fs(x-2, y, z, v, m, s))/(4.d0*dx)       &
          )          &
          ! Term 1
        + Ky * (                                                                            &  
            4.d0/3.d0 * (fs(x, y+1, z, v, m, s) - fs(x, y-1, z, v, m, s))/(2.d0*dy)      &
          - 1.d0/3.d0 * (fs(x, y+2, z, v, m, s) - fs(x, y-2, z, v, m, s))/(4.d0*dy)      &

          ))                                             
        endif

      !&
          ! v_parallel non-linearity
!          - q/(T) * Ez(x,y,z,m,s) * (3./v_scale**2) * ( fs(x, y, z, v, m, s)   &
!          + q/(T) * Ez(x,y,z,m,s) * 1./scale_v (                                         & 
          !            4.d0/3.d0*(fs(x, y, z, v+1, m, s)   &
!          - 1.d0/3.d0*(fs(x, y, z, v+2, m, s)   &
        if(vpar_term) then         
            fts = fts & 
            + q / (mass * scale_v) * Ez(x,y,z,m,s) * (                                         & 
            4.d0/3.d0*(fs(x, y, z, v+1, m, s) - fs(x, y, z, v-1, m, s))/(2.d0*dv)               &
          - 1.d0/3.d0*(fs(x, y, z, v+2, m, s) - fs(x, y, z, v-2, m, s))/(4.d0*dv)      &
          ) 
          endif

      if (ptrap_term) then
            fts = fts + &
            alpha_s / 2. * Mu(m) * B0**2 * epsilon_t * sin(Zv(z))* &
           (4.d0/3.d0*(fs(x, y, z, v+1, m, s) - fs(x, y, z, v-1, m, s))/(2.d0*dv)               &
          - 1.d0/3.d0*(fs(x, y, z, v+2, m, s) - fs(x, y, z, v-2, m, s))/(4.d0*dv)) 
        endif
        
        
        
        !   time-integrate the distribution function     
        if(rk_step == 1) then
            ft(x, y, z, v, m, s) = fts
        else if((rk_step == 2) .or. (rk_step == 3)) then
            ft(x, y, z, v, m, s) = ft(x, y, z, v, m, s) + 2.d0*fts
        else    
            fts = ft(x, y, z, v, m, s) + fts
        end if
        
        fss(x, y, z, v, m, s) = f(x, y, z, v, m, s) + fts*dt
          
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    return
end subroutine morinishi_sa


subroutine Fortran_Vlasov_EM_Linear(f0, f, fs, fss, ft, H, Xi, phi, Ap, dt, rk_step,&
sp, Vel, Mu, B0, beta, eps, max_Xi, eta_hyper_visc) bind(C, name="Fortran_Vlasov_EM_Linear")

      implicit none

      integer :: rk_step
      double precision  ::  q, mass, sub = 1.d0/2.d0, B0, w_n, w_T, scale_v, scale_n, sigma, alpha, &
                            beta, eps, dfs_dv, dH_dz, dfs_dz,T, N, BoBpar
      double precision  ::  Vel(NvGlB:NvGuB),  Mu(NmGlB:NmGuB), &
                            phi(NxLlB-2:NxLuB+2, NyLlB-2:NyLuB+2, NzLlB:NzLuB, NmLlB:NmLuB, NsLlB:NsLuB),                       &
                            Ap (NxLlB-2:NxLuB+2, NyLlB-2:NyLuB+2, NzLlB:NzLuB, NmLlB:NmLuB, NsLlB:NsLuB),                       &
                            Xi (NxLlB-2:NxLuB+2, NyLlB-2:NyLuB+2, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
                            H(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
                            f(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),           &
                            f0(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),          &
                            fs(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),          &
                            fss(NxLlB:NxLuB, NyLlB:NyLuB, NzLlB:NzLuB, NvLlB:NvLuB, NmLlD:NmLuD, NsLlD:NsLuD),         &
                            ft(NxLlD:NxLuD, NyLlD:NyLuD, NzLlD:NzLuD, NvLlD:NvLuD, NmLlD:NmLuD, NsLlD:NsLuD), fts, dt

      !double precision :: Xi_z, Xi_y(NyLlD:NyLuD), Xi_y_B(NyLlD:NyLuD), cs_a=4./3., cs_b=-1./3., cs_c=0., cs_alpha=0.
      !double precision :: Xi_z, Xi_y(NyLlD:NyLuD), Xi_y_B(NyLlD:NyLuD), cs_a, cs_b, cs_c, cs_alpha=1./3.
      double precision :: Xi_z, Xi_y(NyLlD:NyLuD), Xi_y_B(NyLlD:NyLuD), cs_a, cs_b, cs_c, cs_alpha=3./8.
      integer :: v, z, y, x, s, m
    
      double precision :: max_Xi(0:2), eta_hyper_visc, hyper_visc, phi_y, phi_z
      
      type(species2) :: sp(0:8)
      
      cs_a = 1./6. * ( cs_alpha + 9.)
      cs_b = 1./15.* ( 32. * cs_alpha - 9.)
      cs_c = 1./10. *(- 3.* cs_alpha + 1.)
      !cs_c = 0. !1./10. *(- 3.* cs_alpha + 1.)
      max_Xi(:) = 0.
   
     if(do_gyro .eqv. .TRUE.) sub = 3.0d0/2.0d0
         

     do s = NsLlD, NsLuD
       ! charge over mass
       q        = sp(s)%charge
       mass     = sp(s)%mass
       w_n      = sp(s)%w_n
       w_T      = sp(s)%w_T
       T        = sp(s)%T0
       n        = sp(s)%n0
       scale_v  = sp(s)%scale_v
       scale_n  = sp(s)%scale_n
       sigma    = sp(s)%sigma
       alpha    = sp(s)%alpha
    
     !  write(*,*) q,mass, w_n, w_T, T, n, scale_v, scale_n, sigma, alpha, B0
 
        do m = NmLlB, NmLuB
            do v = NvLlB, NvLuB
                do  z = NzLlB, NzLuB
                    do  y = NyLlB-2, NyLuB+2
                        do  x = NxLlB-2, NxLuB+2
                         ! Electro-static 
                            Xi(x,y,z,v,m,s) = phi(x,y,z,m,s)
                         ! Electro-magnetic with Ap
                         !   Xi(x,y,z,v,m,s) = phi(x,y,z,m,s)  - alpha * eps * beta * Vel(v)*Ap(x,y,z,m,s)
                         ! Electro-magnetic with Ap and Bp
                         !   Xi(x,y,z,v,m,s) = phi(x,y,z,m,s)  - alpha * eps * beta * Vel(v)*Ap(x,y,z,m,s) + Mu(m) / q * Bp(x,y,z,m,s)
                        enddo
                    enddo
                    
                    do  y = NyLlB, NyLuB
                        do  x = NxLlB, NxLuB
                            !H(x,y,z,v,m,s)  = fs(x,y,z,v,m,s) + sigma * Xi(x,y,z,v,m,s) * f0(x,y,z,v,m,s)
                            H(x,y,z,v,m,s)  = fs(x,y,z,v,m,s) - sigma * Xi(x,y,z,v,m,s) * f0(x,y,z,v,m,s)
                        enddo
                    enddo
                enddo
            enddo
        enddo
      enddo
       
 dfs_dz = 0.
 hyper_visc = 0.
     do s = NsLlD, NsLuD
       ! charge over mass
       q     = sp(s)%charge
       mass  = sp(s)%mass
       w_n   = sp(s)%w_n
       w_T   = sp(s)%w_T
       T     = sp(s)%T0
       n     = sp(s)%n0
       scale_v    = sp(s)%scale_v
       scale_n    = sp(s)%scale_n

        do m = NmLlD, NmLuD
            do v = NvLlD, NvLuD
                do  z = NzLlD, NzLuD
                        do  x = NxLlD, NxLuD
                    
                        do  y = NyLlD, NyLuD
                            !! calcualte coefficients for compact stencil
                            Xi_y_B(y) = ( cs_a *(Xi(x  , y+1, z, v, m, s) - Xi(x  , y-1, z, v, m, s))/2.  &
                                      +   cs_b *(Xi(x  , y+2, z, v, m, s) - Xi(x  , y-2, z, v, m, s))/4.  &    
                                      +   cs_c *(Xi(x  , y+3, z, v, m, s) - Xi(x  , y-3, z, v, m, s))/6.)/(dy)    
                        enddo
                        
                        call solveCycCCSymTriDiagonal(Xi_y_B, Xi_y, cs_alpha, NyLuD-NyLlD+1)
                        
                        do  y = NyLlD, NyLuD

                        ! Magnetic pre-factor
                        !BoBpar  = 1./(1+beta * sqrt(m*T/2. * j(x,y,z,m,s)/q*B**2*Vel(v)))

                    
                      ! calculate here the derivatives otherwise we would have to allocate whole array, how (in)efficient is this ?
                      ! calculate \frac{\Xe}{x} 
                      ! probably we can improve this by using cycling shift ?
                      ! 6th-order
                      !Xi_y(y)     = (Xi(x  , y+1, z, v, m, s) - Xi(x  , y-1, z, v, m, s))/(2.*dy)    
!                      Xi_y(y)     = (8. *(Xi(x  , y+1, z, v, m, s) - Xi(x  , y-1, z, v, m, s))  &
!                                    -1. *(Xi(x  , y+2, z, v, m, s) - Xi(x  , y-2, z, v, m, s)))/(12.*dy)    
                      
                    !   Xi_y(y)     = (45. *(Xi(x  , y+1, z, v, m, s) - Xi(x  , y-1, z, v, m, s))  &
                    !               -9. *(Xi(x  , y+2, z, v, m, s) - Xi(x  , y-2, z, v, m, s))  &    
                    !               +1. *(Xi(x  , y+3, z, v, m, s) - Xi(x  , y-3, z, v, m, s)))/(60.*dy)    

			if (((w_n + w_T * (((Vel(v)**2.e0+Mu(m)* B0)  - sub)))  * f0(x, y, z, v, m, s )  >= 0.)) then
!                        	Xi_y(y)  = (Xi(x,y-2,z,v,m,s) - 6. * Xi(x,y-1,z,v,m,s) + 3.*Xi(x,y,z,v,m,s) + 2.*Xi(x,y+1,z,v,m,s))/(6.*dy)
!                                Xi_y(y)  = (-Xi(x,y-1,z,v,m,s) + Xi(x,y,z,v,m,s))/(dy)
!                               Xi_y(y)  =  (-2.*Xi(x,y-3,z,v,m,s) + 15. * Xi(x,y-2,z,v,m,s) - 60.*Xi(x,y-1,z,v,m,s) + 20.*Xi(x,y,z,v,m,s) &
!					+30.*Xi(x,y+1,z,v,m,s) - 3. * Xi(x,y+2,z,v,m,s))/(60.*dy)
			else
!                        	Xi_y(y)  = -(Xi(x,y+2,z,v,m,s) - 6. * Xi(x,y+1,z,v,m,s) + 3.*Xi(x,y,z,v,m,s) + 2.*Xi(x,y-1,z,v,m,s))/(6.*dy)
!                                Xi_y(y)  = -(-Xi(x,y+1,z,v,m,s) + Xi(x,y,z,v,m,s))/(dy)
!                               Xi_y(y)  = -(-2.*Xi(x,y+3,z,v,m,s) + 15. * Xi(x,y+2,z,v,m,s) - 60.*Xi(x,y+1,z,v,m,s) + 20.*Xi(x,y,z,v,m,s) &
!					+30.*Xi(x,y-1,z,v,m,s) - 3. * Xi(x,y-2,z,v,m,s))/(60.*dy)
			endif

                       Xi_z     = (8. *(Xi(x, y, z+1, v, m, s) - Xi(x, y, z-1, v, m, s)) &
                                  -1. *(Xi(x, y, z+2, v, m, s) - Xi(x, y, z-2, v, m, s)))/(12.*dz) 
                      
                       dfs_dv   = (8. *(fs(x, y, z, v+1, m, s) - fs(x, y, z, v-1, m, s))               &
                                  -1. *(fs(x, y, z, v+2, m, s) - fs(x, y, z, v-2, m, s)))/(12.*dv)

                       dH_dz   = (8.*(H (x, y, z+1, v, m, s) - H (x, y, z-1, v, m, s))               &
                                 -1.*(H (x, y, z+2, v, m, s) - H (x, y, z-2, v, m, s)))/(12.*dz)

                       !dfs_dz  = (fs(x,y,z+1, v,m,s) - fs(x,y,z-1,v,m,s))/(2.*dz)
                       ! 3rd order upwind (6*dx)^(-1) [ 1 -6 3 2 0 ]
                       !if (Vel(v) >= 0.) then
                       !dfs_dz  = (-fs(x,y,z+2, v,m,s) + 6. * fs(x,y,z+1,v,m,s) - 3.*fs(x,y,z,v,m,s) - 2.*fs(x,y,z-1,v,m,s))/(6.*dz)
                       !else
!                            dfs_dz  = (fs(x,y,z-2, v,m,s) - 6. * fs(x,y,z-1,v,m,s) + 3.*fs(x,y,z,v,m,s) + 2.*fs(x,y,z+1,v,m,s))/(6.*dz)
                       !endif
        
                        !hyper_visc_z = eta_hyper_visc * (-fs(x,y,z-2,v,m,s)+4.* fs(x,y,z-1,v,m,s)-6.*fs(x,y,z,v,m,s)+4.*fs(x,y,z+1,v,m,s)-fs(x,y,z+2,v,m,s))/(dz**4)
                        hyper_visc = 0. !eta_hyper_visc * (fs(x,y-2,z,v,m,s)-4.* fs(x,y-1,z,v,m,s)+6.*fs(x,y,z,v,m,s)-4.*fs(x,y+1,z,v,m,s)+fs(x,y+2,z,v,m,s))/(dy**4)

                        ! add hyper-viscosity to stabilize high k_parallel modes 
                        ! Ref : Pueschel, Dannert, Jenko (Comp. Phys. Comm.) 2010 : On the role of numerical Dissipation in gyrokinetic Vlasov simulations of plasma microturbulence
              !          hyper_visc = 0. !eta_hyper_visc * (fs(x,y,z-2,v,m,s)-4.* fs(x,y,z-1,v,m,s)+6.*fs(x,y,z,v,m,s)-4.*fs(x,y,z+1,v,m,s)+fs(x,y,z+2,v,m,s))/16.


                        if(rk_step == 4) then 
                            max_Xi(1) = max(max_Xi(1), abs(Xi_y(y)))
                            max_Xi(2) = max(max_Xi(2), abs(Xi_z))
                        endif
                      phi_y     = (8. *(phi(x  , y+1, z, m, s) - phi(x  , y-1, z, m , s))  &
                                    -1. *(phi(x  , y+2, z, m, s) - phi(x  , y-2, z, m , s)))/(12.*dy)    
                       phi_z     = (8. *(phi(x, y, z+1, m , s) - phi(x, y, z-1, m , s)) &
                                  -1. *(phi(x, y, z+2, m, s) - phi(x, y, z-2, m, s)))/(12.*dz) 
      ! calculate the time derivatve      
          
           fts = - 4.d0/3.d0*Vel(v)*scale_v*(fs(x, y, z+1, v, m, s)        &
          - fs(x, y, z-1, v, m, s))/(2.d0*dz)               &
          + 1.d0/3.d0*Vel(v)*scale_v*(fs(x, y, z+2, v, m, s)        &
          - fs(x, y, z-2, v, m, s))/(4.d0*dz)               &
           + phi_y* (w_n + w_T * (((Vel(v)**2.e0+Mu(m)* B0)  - sub)))  * f0(x, y, z, v, m, s ) &
           + q/T*Vel(v)*scale_v  * phi_z* f0(x, y, z, v, m, s)   
                        
                        ! Finally the Vlasov equation calculate the time derivatve      
       !                 fts =  -alpha * Vel(v)* dH_dz   &
       !                       + Xi_y(y) * (w_n + w_T * (Vel(v)**2.e0+Mu(m)* B0  - sub))  * f0(x, y, z, v, m, s )   !&
                    !    fts = - alpha * Vel(v)* dH_dz   &
                    !          - Xi_y(y) * (w_n + w_T * (Vel(v)**2.e0+Mu(m)* B0  - sub))  * f0(x, y, z, v, m, s )   !&
                        
                              !          + sigma/mass/scale_v**2 * Xi_z* dfs_dv
                            !- Xi_y(y) * (w_n + w_T * (Vel(v)**2.e0+Mu(m)* B0  - sub)  * f0(x, y, z, v, m, s ) + hyper_visc  !&

        !   time-integrate the distribution function     
        if(rk_step == 1) then
            ft(x, y, z, v, m, s) = fts
        else if((rk_step == 2) .or. (rk_step == 3)) then
            ft(x, y, z, v, m, s) = ft(x, y, z, v, m, s) + 2.d0*fts
        else    
            fts = ft(x, y, z, v, m, s) + fts
        end if
        
        fss(x, y, z, v, m, s) = f(x, y, z, v, m, s) + fts*dt
          
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    return
end subroutine Fortran_Vlasov_EM_Linear


end module Morinishi
