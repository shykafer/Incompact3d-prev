!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module time_integrators

  implicit none

  private
  public :: int_time, set_time_int_coefficient

contains
  !############################################################################
  subroutine set_time_int_coefficient()
    use decomp_2d, only : mytype
    use param, only : adt,bdt,cdt,ddt,gdt,dt
    use param, only : ntime,nrhotime,itimescheme,iadvance_time
    use param, only : zero, one
 
    implicit none 

    !
    adt(:)=zero ; bdt(:)=zero ; cdt(:)=zero ; gdt(:)=zero
    if (itimescheme.eq.1) then ! Euler
       iadvance_time=1
       adt(1)=1.0_mytype*dt
       bdt(1)=0.0_mytype*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)
 
       ntime = 1
       nrhotime = 2
    elseif (itimescheme.eq.2) then ! AB2
       iadvance_time=1
       adt(1)=1.5_mytype*dt
       bdt(1)=-0.5_mytype*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)
 
       ntime = 2
       nrhotime = 3
    elseif (itimescheme.eq.3) then ! AB3
       iadvance_time=1
       adt(1)= (23._mytype/12._mytype)*dt
       bdt(1)=-(16._mytype/12._mytype)*dt
       cdt(1)= ( 5._mytype/12._mytype)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)
 
       ntime = 3
       nrhotime = 4
    elseif(itimescheme==4) then  ! AB4
       iadvance_time=1
       adt(1)=(55.0_mytype/24.0_mytype)*dt
       bdt(1)=-(59.0_mytype/24.0_mytype)*dt
       cdt(1)=(37.0_mytype/24.0_mytype)*dt
       ddt(1)=-(9.0_mytype/24.0_mytype)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)+ddt(1)
       gdt(3)=gdt(1)
 
       ntime = 4
       nrhotime = 5
    elseif(itimescheme.eq.5) then !RK3
       iadvance_time=3
       adt(1)=(8._mytype/15._mytype)*dt
       bdt(1)=0._mytype
       gdt(1)=adt(1)
       adt(2)=(5._mytype/12._mytype)*dt
       bdt(2)=(-17._mytype/60._mytype)*dt
       gdt(2)=adt(2)+bdt(2)
       adt(3)=(3._mytype/4._mytype)*dt
       bdt(3)=(-5._mytype/12._mytype)*dt
       gdt(3)=adt(3)+bdt(3)
 
       ntime = 2
       nrhotime = 3
    elseif(itimescheme.eq.6) then !RK4 Carpenter and Kennedy
       iadvance_time=5
       adt(1)=0.0_mytype
       adt(2)=-0.4178904745_mytype
       adt(3)=-1.192151694643_mytype
       adt(4)=-1.697784692471_mytype
       adt(5)=-1.514183444257_mytype
       bdt(1)=0.1496590219993_mytype
       bdt(2)=0.3792103129999_mytype
       bdt(3)=0.8229550293869_mytype
       bdt(4)=0.6994504559488_mytype
       bdt(5)=0.1530572479681_mytype
       gdt(1)=0.1496590219993_mytype*dt
       gdt(2)=0.220741935365_mytype*dt
       gdt(3)=0.25185480577_mytype*dt
       gdt(4)=0.33602636754_mytype*dt
       gdt(5)=0.041717869325_mytype*dt
 
       ntime = 2
       nrhotime = 5 ! (A guess)
 
    elseif(itimescheme.eq.7) then !Semi-implicit
       iadvance_time=1
       adt(1)= (23./12.)*dt
       bdt(1)=-(16./12.)*dt
       cdt(1)= ( 5./12.)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)
 
       ntime = 3
       nrhotime = 4
    endif
 
  end subroutine set_time_int_coefficient
  !############################################################################

  subroutine intt(var1,dvar1,forcing1)

    USE param
    USE variables
    USE decomp_2d
    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), optional :: forcing1

#ifdef DEBG
    if (nrank .eq. 0) print *,'# intt start'
#endif

    if (itimescheme.eq.1) then
       !>>> Euler
       var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
    elseif(itimescheme.eq.2) then
       !>>> Adam-Bashforth second order (AB2)

       ! Do first time step with Euler
       if(itime.eq.1.and.irestart.eq.0) then
          var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
       else
          var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
       endif
       dvar1(:,:,:,2)=dvar1(:,:,:,1)
    elseif(itimescheme.eq.3) then
       !>>> Adams-Bashforth third order (AB3)

       ! Do first time step with Euler
       if(itime.eq.1.and.irestart.eq.0) then
          var1(:,:,:)=dt*dvar1(:,:,:,1)+var1(:,:,:)
       elseif(itime.eq.2.and.irestart.eq.0) then
          ! Do second time step with AB2
          var1(:,:,:)=onepfive*dt*dvar1(:,:,:,1)-half*dt*dvar1(:,:,:,2)+var1(:,:,:)
          dvar1(:,:,:,3)=dvar1(:,:,:,2)
       else
          ! Finally using AB3
          var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+cdt(itr)*dvar1(:,:,:,3)+var1(:,:,:)
          dvar1(:,:,:,3)=dvar1(:,:,:,2)
       endif
       dvar1(:,:,:,2)=dvar1(:,:,:,1)
    elseif(itimescheme.eq.4) then
       !>>> Adams-Bashforth fourth order (AB4)

       if (nrank.eq.0) then
          print *, "AB4 not implemented!"
          STOP
       endif

       !if (itime.eq.1.and.ilit.eq.0) then
       !var(:,:,:)=gdt(itr)*hx(:,:,:)+var(:,:,:)
       !uy(:,:,:)=gdt(itr)*hy(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=gdt(itr)*hz(:,:,:)+uz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !elseif (itime.eq.2.and.ilit.eq.0) then
       !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+var(:,:,:)
       !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+uz(:,:,:)
       !gox(:,:,:)=gx(:,:,:)
       !goy(:,:,:)=gy(:,:,:)
       !goz(:,:,:)=gz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !elseif (itime.eq.3.and.ilit.eq.0) then
       !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+var(:,:,:)
       !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+uz(:,:,:)
       !gox(:,:,:)=gx(:,:,:)
       !goy(:,:,:)=gy(:,:,:)
       !goz(:,:,:)=gz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !else
       !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+ddt(itr)*gax(:,:,:)+var(:,:,:)
       !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+ddt(itr)*gay(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+ddt(itr)*gaz(:,:,:)+uz(:,:,:)
       !gax(:,:,:)=gox(:,:,:)
       !gay(:,:,:)=goy(:,:,:)
       !gaz(:,:,:)=goz(:,:,:)
       !gox(:,:,:)=gx(:,:,:)
       !goy(:,:,:)=gy(:,:,:)
       !goz(:,:,:)=gz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !endif
       !>>> Runge-Kutta (low storage) RK3
    elseif(itimescheme.eq.5) then
       if(itr.eq.1) then
          var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
       else
          var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
       endif
       dvar1(:,:,:,2)=dvar1(:,:,:,1)
       !>>> Runge-Kutta (low storage) RK4
    elseif(itimescheme.eq.6) then

       if (nrank.eq.0) then
          print *, "RK4 not implemented!"
          STOP
       endif
       !>>> Semi-implicit
    elseif(itimescheme.eq.7) then

       call inttimp(var1,dvar1,forcing1)

    else

       if (nrank.eq.0) then
          print *, "Unrecognised itimescheme: ", itimescheme
          STOP
       endif

    endif

#ifdef DEBG
    if (nrank .eq. 0) print *,'# intt done'
#endif

    return

  end subroutine intt

  SUBROUTINE int_time(rho1, ux1, uy1, uz1, phi1, drho1, dux1, duy1, duz1, dphi1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : zero, one
    USE param, ONLY : ntime, nrhotime, ilmn, iscalar, ilmn_solve_temp,itimescheme
    USE param, ONLY : primary_species, massfrac
    use param, only : scalar_lbound, scalar_ubound
    USE variables, ONLY : numscalar,nu0nu
    USE var, ONLY : ta1, tb1

    IMPLICIT NONE

    !! INPUT/OUTPUT
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1, dux1, duy1, duz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! OUTPUT
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

    !! LOCAL
    INTEGER :: is, i, j, k

    CALL int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

    IF (ilmn) THEN
       IF (ilmn_solve_temp) THEN
          CALL int_time_temperature(rho1, drho1, dphi1, phi1)
       ELSE
          CALL int_time_continuity(rho1, drho1)
       ENDIF
    ENDIF

    IF (iscalar.NE.0) THEN
       IF (ilmn.and.ilmn_solve_temp) THEN
          !! Compute temperature
          call calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))
       ENDIF

       DO is = 1, numscalar
          IF (is.NE.primary_species) THEN
             IF (itimescheme.ne.7) THEN
                CALL intt(phi1(:,:,:,is), dphi1(:,:,:,:,is))
             ELSE
            !!TO BE DONE: when sc is not one
                CALL scalarimp(ux1,uy1,uz1,phi1(:,:,:,is),dphi1(:,:,:,:,is),is)
             ENDIF

             DO k = 1, xsize(3)
                DO j = 1, xsize(2)
                   DO i = 1, xsize(1)
                      phi1(i,j,k,is) = max(phi1(i,j,k,is),scalar_lbound(is))
                      phi1(i,j,k,is) = min(phi1(i,j,k,is),scalar_ubound(is))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       IF (primary_species.GE.1) THEN
          phi1(:,:,:,primary_species) = one
          DO is = 1, numscalar
             IF ((is.NE.primary_species).AND.massfrac(is)) THEN
                phi1(:,:,:,primary_species) = phi1(:,:,:,primary_species) - phi1(:,:,:,is)
             ENDIF
          ENDDO

          DO k = 1, xsize(3)
             DO j = 1, xsize(2)
                DO i = 1, xsize(1)
                   phi1(i,j,k,primary_species) = max(phi1(i,j,k,primary_species),zero)
                   phi1(i,j,k,primary_species) = min(phi1(i,j,k,primary_species),one)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       IF (ilmn.and.ilmn_solve_temp) THEN
          !! Compute rho
          call calc_rho_eos(rho1(:,:,:,1), ta1, phi1, tb1, xsize(1), xsize(2), xsize(3))
       ENDIF
    ENDIF

  ENDSUBROUTINE int_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_momentum
  !! DESCRIPTION: Integrates the momentum equations in time by calling time
  !!              integrator.
  !!      INPUTS: dux1, duy1, duz1 - the RHS(s) of the momentum equations
  !!     OUTPUTS: ux1,   uy1,  uz1 - the intermediate momentum state.
  !!       NOTES: This is integrating the MOMENTUM in time (!= velocity)
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

    USE param
    USE variables
    USE var, ONLY: px1, py1, pz1
    USE decomp_2d

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1, duy1, duz1

    if (itimescheme.eq.7) then
       call intt(ux1, dux1, px1)
       call intt(uy1, duy1, py1)
       call intt(uz1, duz1, pz1)
    else
       call intt(ux1, dux1)
       call intt(uy1, duy1)
       call intt(uz1, duz1)
    endif

  endsubroutine int_time_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_continuity
  !! DESCRIPTION: Integrates the continuity (aka density transport) equation in
  !!              time
  !!      INPUTS: drho1 - the RHS(s) of the continuity equation.
  !!     OUTPUTS:  rho1 - the density at new time.
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_continuity(rho1, drho1)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    integer :: it, i, j, k
    real(mytype) :: rhomin, rhomax

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! First, update old density / store old transients depending on scheme
    if (itimescheme.lt.5) then
       !! Euler/AB - Store old density values
       do it = nrhotime, 2, -1
          rho1(:,:,:,it) = rho1(:,:,:,it-1)
       enddo
    elseif (itimescheme.eq.5) then
       !! RK3 - Stores old transients
       if (itr.eq.1) then
          do it = nrhotime, 2, -1
             rho1(:,:,:,it) = rho1(:,:,:,it-1)
          enddo
          rho1(:,:,:,2) = drho1(:,:,:,1)
       endif
    else
       if (nrank.eq.0) then
          print *, "int_time_continuity not implemented for itimescheme", itimescheme
          stop
       endif
    endif

    !! Now we can update current density
    call intt(rho1(:,:,:,1), drho1)

    !! Enforce boundedness on density
    if (ilmn_bound) then
       rhomin = min(dens1, dens2)
       rhomax = max(dens1, dens2)
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                rho1(i, j, k, 1) = max(rho1(i, j, k, 1), rhomin)
                rho1(i, j, k, 1) = min(rho1(i, j, k, 1), rhomax)
             enddo
          enddo
       enddo
    endif

  endsubroutine int_time_continuity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_temperature
  !! DESCRIPTION: Integrates the temperature equation in time
  !!      INPUTS: drho1 - the RHS(s) of the temperature equation.
  !!     OUTPUTS:  rho1 - the density at new time.
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_temperature(rho1, drho1, dphi1, phi1)

    USE param
    USE variables
    USE decomp_2d

    USE navier, ONLY : lmn_t_to_rho_trans
    USE var, ONLY : tc1, tb1

    implicit none

    integer :: it, i, j, k

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! First, update old density / store old transients depending on scheme
    if (itimescheme.lt.5) then
       !! Euler/AB - Store old density values
       do it = nrhotime, 2, -1
          rho1(:,:,:,it) = rho1(:,:,:,it-1)
       enddo
    elseif (itimescheme.eq.5) then
       !! RK3 - Stores old transients
       if (itr.eq.1) then
          do it = nrhotime, 2, -1
             rho1(:,:,:,it) = rho1(:,:,:,it-1)
          enddo

          !! Convert temperature transient to density transient and store it.
          call lmn_t_to_rho_trans(rho1(:,:,:,2), drho1(:,:,:,1), rho1(:,:,:,1), dphi1, phi1)
       endif
    else
       if (nrank.eq.0) then
          print *, "int_time_continuity not implemented for itimescheme", itimescheme
          stop
       endif
    endif

    !!-------------------------------------------------------------------
    !! XXX We are integrating the temperature equation - get temperature
    !!-------------------------------------------------------------------
    call calc_temp_eos(tc1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

    !! Now we can update current temperature
    call intt(tc1, drho1)

    !! Temperature >= 0
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             tc1(i,j,k) = max(tc1(i,j,k), zero)
          enddo
       enddo
    enddo

    !!-------------------------------------------------------------------
    !! XXX We are integrating the temperature equation - get back to rho
    !!-------------------------------------------------------------------
    call calc_rho_eos(rho1(:,:,:,1), tc1, phi1, tb1, xsize(1), xsize(2), xsize(3))

  endsubroutine int_time_temperature

end module time_integrators
!##################################################################
!##################################################################
module time_step

   implicit none
 
   private
   public :: time_stepping_init, time_stepping

   contains
   !############################################################################
   subroutine time_stepping_init()
      use param, only : cfl_diff_sum, cfl_diff_x, cfl_diff_y, cfl_diff_z
      use param, only : dt, dtmax, diff_crit, dtstep
      use decomp_2d, only : mytype, nrank
      implicit none
      
      real(mytype) :: diffmax

      ! Check diffusion number -> maximum time step
      diffmax = maxval((/cfl_diff_x, cfl_diff_y, cfl_diff_z /))
      dtmax = floor((diff_crit/diffmax * dt)/dtstep)*dtstep  ! Control round off errors by dtstep

      if (nrank==0) write(*,"(' dtmax             : ',F17.10)") dtmax
      !write(*,"(' rank',I4,' dtmax             : ',F17.10)") nrank,dtmax

   end subroutine time_stepping_init
   !##################################################################
   !############################################################################
   !!
   !!  SUBROUTINE: time_stepping
   !! DESCRIPTION: Sets new time step for fixed or dynamic time step.
   !!      AUTHOR: Kay Schäfer
   !!
   !############################################################################
   subroutine time_stepping(itime,ux1,uy1,uz1,ep1)
      use tools, only : compute_cfl
      use time_integrators, only : set_time_int_coefficient

      use param, only : itime0, t0, t, dt, ifirst, idyndt, fdyndt, iibm
      use param, only : cfl_crit, diff_crit, dtmax, dtstep, cflx_max, cfly_max, cflz_max
      use param, only : zero, one
      use var, only : xsize
      use decomp_2d, only : mytype, nrank
      use ibm, only : body
      
      implicit none 

      !! INPUTS
      integer :: itime
      real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1

      real(mytype) :: cflmax
      
      if (idyndt.eq.zero) then 
         t=t0 + (itime0 + itime + 1 - ifirst)*dt
      else if (idyndt.eq.one) then 
         if ((mod(itime, fdyndt).eq.0).or.(itime.eq.ifirst)) then
            ! In case of IBM set values in body to zero
            if (iibm.ne.0) call body(ux1,uy1,uz1,ep1)

            ! compute CFL numbers
            call compute_cfl(ux1,uy1,uz1,1)
            cflmax = maxval((/cflx_max, cfly_max, cflz_max /))

            dt = floor((cfl_crit/cflmax * dt)/dtstep)*dtstep  ! Control round off errors by dtstep
            !write(*,"(' rank',I4,' dt before        : ',F17.10)") nrank,dt
            if (dt.gt.dtmax) dt = dtmax 
            !write(*,"(' rank',I4,' dt after         : ',F17.10)") nrank,dt
            
            ! update time integration coefficients
            call set_time_int_coefficient()
         end if 
         ! update time 
         t = t + dt
      end if
   
   endsubroutine time_stepping
   !##################################################################
end module time_step
