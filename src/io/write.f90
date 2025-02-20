!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module write_sub
  use math_constants,only : zi,pi
  implicit none

contains

!===================================================================================================================================

  !! export SYSNAME_k.data file
  subroutine write_k_data(system,stencil)
    use structures
    use salmon_global, only: base_directory,sysname,yn_periodic
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root,comm_sync_all
    use filesystem, only: open_filehandle
    implicit none
    type(s_dft_system) ,intent(in) :: system
    type(s_stencil),intent(in) :: stencil
    !
    integer :: fh_k
    integer :: ik,NK
    character(100) :: file_k_data
    
    if(yn_periodic=='n') then
      return
    end if

    NK = system%nk
    file_k_data = trim(base_directory)//trim(sysname)//'_k.data'

    if (comm_is_root(nproc_id_global)) then
      fh_k = open_filehandle(file_k_data, status="replace")
      if(stencil%if_orthogonal) then
        write(fh_k, '("#",1X,A)') "k-point distribution"
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "Reduced coordinate of k-points"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "none", &
          & 3, "ky", "none", &
          & 4, "kz", "none", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,E23.15E3))') &
            & ik, &
            & system%vec_k(1,ik) / system%primitive_b(1,1), &
            & system%vec_k(2,ik) / system%primitive_b(2,2), &
            & system%vec_k(3,ik) / system%primitive_b(3,3), &
            & system%wtk(ik)
        end do !ik
        write(fh_k, '("#",1X,A)') "coefficients (2*pi/a [a.u.]) in kx, ky, kz"
        write(fh_k, '(3E23.15E3)') system%primitive_b(1,1), system%primitive_b(2,2), system%primitive_b(3,3) 

      else
        write(fh_k, '("#",1X,A)') "k-point distribution (nonorthogonal coordinate)"
        write(fh_k, '("#",1X,A,":",1X,A)') "brl", "reciprocal primitive vectors"
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "k-vectors"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        write(fh_k, '("#",A,"[",A,"]")') "brl(1:3,1:3)", "a.u."
        write(fh_k, '(9(1X,E23.15E3))') system%primitive_b(1:3,1:3)
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "a.u.", &
          & 3, "ky", "a.u.", &
          & 4, "kz", "a.u.", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,E23.15E3))') &
            & ik, &
            & system%vec_k(1,ik), &
            & system%vec_k(2,ik), &
            & system%vec_k(3,ik), &
            & system%wtk(ik)
        end do !ik
      end if
      close(fh_k)
    end if
    call comm_sync_all
  end subroutine write_k_data

!===================================================================================================================================

  !! export SYSNAME_tm.data, SYSNAME_sigma.data, SYSNAME_epsilon.data file
  subroutine write_tm_data(tpsi,system,info,mg,stencil,srg,ppg,energy,yn_ikr)
    use structures
    use stencil_sub
    use sendrecv_grid
    use salmon_global, only: yn_out_tm,yn_out_gs_sgm_eps, &
                       out_gs_sgm_eps_mu_nu, out_gs_sgm_eps_width, &
                       base_directory,sysname, de,nenergy,nelec,xc
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root,comm_summation,comm_sync_all
    use filesystem, only: open_filehandle
    use inputoutput, only: t_unit_energy
    implicit none
    type(s_dft_system) ,intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital)       :: tpsi
    type(s_orbital)       :: tpsi_backup
    type(s_dft_energy) :: energy
    character(1),intent(in) :: yn_ikr
    !
    logical :: flag_print_tm, flag_print_eps
    integer :: fh_tm, narray
    integer :: i,j,ik,ib,ib1,ib2,ilma,nlma,ia,ix,iy,iz,NB,NK,im,ispin
    integer :: ik_s,ik_e,io_s,io_e,is(3),ie(3)
    real(8) :: x,y,z
    complex(8),allocatable :: upu(:,:,:,:),upu_l(:,:,:,:)
    complex(8),allocatable :: upu_all(:,:,:,:),upu_all_l(:,:,:,:)
    complex(8),allocatable :: gtpsi(:,:,:,:),gtpsi_l(:,:,:,:)
    complex(8),allocatable :: uVpsi(:),uVpsi_l(:),uVrpsi(:,:),uVrpsi_l(:,:)
    complex(8),allocatable :: u_rVnl_Vnlr_u(:,:,:,:)
    complex(8),allocatable :: u_rVnl_Vnlr_u_all(:,:,:,:),u_rVnl_Vnlr_u_all_l(:,:,:,:)
    complex(8) :: u_rVnl_u(3),u_Vnlr_u(3),veik
    complex(8) :: wrk(3)
    character(100) :: file_tm_data
    !(for printing dielectric function)
    integer :: mu,nu,iw,nomega
    real(8) :: w, omega_max, domega, delta, delta_munu, n_e, V, deigen
    complex(8),allocatable :: matrix_vij(:,:,:,:), matrix_vji(:,:,:,:)
    complex(8),allocatable :: mat_l(:,:), mat(:,:)
    complex(8) :: sigma,eps,sigma_intra,sigma_inter,eps_intra,eps_inter
    complex(8) :: sigma_l,sigma_intra_l
    character(256) :: filename
    integer,parameter :: fh_s=333,fh_e=777
    complex(8),allocatable :: expikr(:,:,:,:)

    flag_print_tm = .false.
    flag_print_eps= .false.
    if(yn_out_tm=='y') flag_print_tm = .true.
    if(yn_out_gs_sgm_eps=='y') flag_print_eps= .true.

    !For yn_out_gs_sgm_eps=='y'
    !following input parameters are used:
    !   out_gs_sgm_eps_mu_nu(1:2) !=mu,nu (x,y,z)
    !   out_gs_sgm_eps_width      != delta (width)
    !   nenergy         != # of eneryg(omega) point
    !   de              != dw




    if(info%im_s/=1 .or. info%im_e/=1) then!??????
      write(*,*) "error @ write_tm_data: im/=1"
      return
    endif
    if(system%Nspin/=1) then!??????
      write(*,*) "error @ write_tm_data: nspin/=1"
      return
    endif
    !if(info%io_s/=1 .or. info%io_e/=system%no) then!??????
    !  if (comm_is_root(nproc_id_global)) then
    !    write(*,*) "error @ write_tm_data: do not use orbital parallelization"
    !    write(*,*) "Only <u|p|u> is printed (pseudopotential terms are not printed)"
    !  endif
    ! !return
    !endif
    if(.not. allocated(tpsi%zwf)) then!??????
      write(*,*) "error @ write_tm_data: do not use real wavefunction (iperiodic=0)"
      return
    endif

    if (comm_is_root(nproc_id_global)) &
    write(*,*) "  calculating transition moment ....."

    im = 1
    ispin = 1

    NB = system%no
    NK = system%nk

    ik_s = info%ik_s
    ik_e = info%ik_e

    io_s = info%io_s
    io_e = info%io_e

    is = mg%is
    ie = mg%ie

    Nlma = ppg%Nlma

    if(yn_ikr=='y')then
       if(.not. allocated(tpsi_backup%zwf)) then
         call allocate_orbital_complex(system%nspin,mg,info,tpsi_backup)
       end if
       tpsi_backup%zwf(:,:,:,:,:,:,:) = tpsi%zwf(:,:,:,:,:,:,:)
   
       allocate(expikr(is(1):ie(1),is(2):ie(2),is(3):ie(3),ik_s:ik_e))
       do ik=ik_s,ik_e
         do iz=is(3),ie(3)
         do iy=is(2),ie(2)
         do ix=is(1),ie(1)
           expikr(ix,iy,iz,ik) = exp(zi*(system%vec_k(1,ik)*dble(ix-1)*system%Hgs(1) &
                                        +system%vec_k(2,ik)*dble(iy-1)*system%Hgs(2) &
                                        +system%vec_k(3,ik)*dble(iz-1)*system%Hgs(3)))
         end do
         end do
         end do
       end do
   
       do ik=ik_s,ik_e
         do ib=io_s,io_e
           do iz=is(3),ie(3)
           do iy=is(2),ie(2)
           do ix=is(1),ie(1)
             tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) = expikr(ix,iy,iz,ik) * tpsi%zwf(ix,iy,iz,ispin,ib,ik,im)
           end do
           end do
           end do
         end do
       end do
       deallocate(expikr)
    end if

    !calculate <u_nk|p_j|u_mk>  (j=x,y,z)

    allocate( upu_l(3,io_s:io_e,NB,ik_s:ik_e) )
    !$omp parallel do private(ik,ib1,ib2) collapse(3)
    do ik =ik_s,ik_e
    do ib1=io_s,io_e
    do ib2=1,NB
       upu_l(1:3,ib1,ib2,ik) = 0d0
    enddo
    enddo
    enddo
    !$omp end parallel do

    allocate(gtpsi_l(3,mg%is_array(1):mg%ie_array(1) &
                      ,mg%is_array(2):mg%ie_array(2) &
                      ,mg%is_array(3):mg%ie_array(3)))
    allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
                    ,mg%is_array(2):mg%ie_array(2) &
                    ,mg%is_array(3):mg%ie_array(3)))
    narray= 3 * ( mg%ie_array(1) - mg%is_array(1) + 1 ) &
              * ( mg%ie_array(2) - mg%is_array(2) + 1 ) &
              * ( mg%ie_array(3) - mg%is_array(3) + 1 )

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if

    do ik=ik_s,ik_e
    do ib2=1,NB

      ! gtpsi = (nabla) psi
      !$omp workshare
      gtpsi_l(:,:,:,:) = 0d0
      !$omp end workshare
      if(ib2.ge.io_s .and. ib2.le.io_e) &
        call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,ib2,ik,im),gtpsi_l,mg%is_array,mg%ie_array,mg%is,mg%ie &
            ,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
      call comm_summation(gtpsi_l,gtpsi,narray,info%icomm_o)

      !$omp parallel do private(ib1,i,wrk)
      do ib1=io_s,io_e
        do i=1,3
          wrk(i) = sum(conjg(tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,ib1,ik,im)) &
                            * gtpsi(i,is(1):ie(1),is(2):ie(2),is(3):ie(3)) )
        end do
        upu_l(:,ib1,ib2,ik) = - zI * wrk * system%Hvol
      end do
      !$omp end parallel do
    end do
    end do
    deallocate(gtpsi_l,gtpsi)


    allocate( upu(3,io_s:io_e,NB,ik_s:ik_e) )
    call comm_summation(upu_l,upu,3*(io_e-io_s+1)*NB*(ik_e-ik_s+1),info%icomm_r)
!    call comm_summation(upu_l,upu,3*NB*NB*(ik_e-ik_s+1),info%icomm_ro)
!    call comm_summation(upu_l,upu,3*NB*NB*NK,info%icomm_rko)
    deallocate(upu_l)


    if(flag_print_tm) then
       allocate( upu_all_l(3,NB,NB,NK), upu_all(3,NB,NB,NK) )
       !$omp parallel do private(ik,ib1,ib2) collapse(3)
       do ik =1,NK
       do ib1=1,NB
       do ib2=1,NB
          if( (ik .ge.ik_s .and. ik .le.ik_e) .and. &
              (ib1.ge.io_s .and. ib1.le.io_e) ) then
             upu_all_l(1:3,ib1,ib2,ik) = upu(1:3,ib1,ib2,ik)
          else
             upu_all_l(1:3,ib1,ib2,ik) = 0d0
          endif
       enddo
       enddo
       enddo
       !$omp end parallel do
       call comm_summation(upu_all_l,upu_all,3*NB*NB*NK,info%icomm_ko)
!       call comm_summation(upu_all_l,upu_all,3*NB*NB*NK,info%icomm_k)
       deallocate(upu_all_l)
    endif


    !calculate -i* <u_mk|[r_j,V_nl]|u_nk>  (j=x,y,z)

    allocate( uVpsi(NB)   ,uVpsi_l(NB) )
    allocate( uVrpsi(3,NB),uVrpsi_l(3,NB) )
    allocate( u_rVnl_Vnlr_u(3,io_s:io_e,NB,ik_s:ik_e) )

    !$omp workshare
    u_rVnl_Vnlr_u(:,:,:,:) = 0d0
    !$omp end workshare

    do ik=ik_s,ik_e
    do ilma=1,Nlma
       ia=ppg%ia_tbl(ilma)
       uVpsi_l  = 0d0
       uVrpsi_l = 0d0
       !$omp parallel do private(j,x,y,z,ix,iy,iz,veik,ib) reduction(+:uVpsi_l,uVrpsi_l)
       do j=1,ppg%Mps(ia)
          x = ppg%Rxyz(1,j,ia)
          y = ppg%Rxyz(2,j,ia)
          z = ppg%Rxyz(3,j,ia)
          ix = ppg%Jxyz(1,j,ia)
          iy = ppg%Jxyz(2,j,ia)
          iz = ppg%Jxyz(3,j,ia)
          veik = conjg(ppg%zekr_uV(j,ilma,ik))
         !do ib=1,NB
          do ib=io_s,io_e
             uVpsi_l(ib)    = uVpsi_l(ib)    + veik*    tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik|u>
             uVrpsi_l(1,ib) = uVrpsi_l(1,ib) + veik* x *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*x|u>
             uVrpsi_l(2,ib) = uVrpsi_l(2,ib) + veik* y *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*y|u>
             uVrpsi_l(3,ib) = uVrpsi_l(3,ib) + veik* z *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*z|u>
          enddo
       end do
       !$omp end parallel do
       uVpsi_l  = uVpsi_l * ppg%rinv_uvu(ilma)
       call comm_summation(uVpsi_l ,uVpsi ,  NB,info%icomm_ro)
       call comm_summation(uVrpsi_l,uVrpsi,3*NB,info%icomm_ro)
!       call comm_summation(uVpsi_l ,uVpsi ,  NB,info%icomm_r)
!       call comm_summation(uVrpsi_l,uVrpsi,3*NB,info%icomm_r)

       !$omp parallel
       !$omp do private(ib1,ib2,u_rVnl_u,u_Vnlr_u) collapse(2)
       do ib1=io_s,io_e
       do ib2=1,NB
          !<u|e^{-ik}*r|v><v|e^{ik}|u>
          u_rVnl_u(1:3) = conjg(uVrpsi(1:3,ib1))*uVpsi(ib2)
        
          !<u|e^{-ik}|v><v|e^{ik}*r|u>
          u_Vnlr_u(1:3) = conjg(uVpsi(ib1))*uVrpsi(1:3,ib2)
        
          u_rVnl_Vnlr_u(:,ib1,ib2,ik) = u_rVnl_Vnlr_u(:,ib1,ib2,ik)  &
          &                           - zi * ( u_rVnl_u(:) - u_Vnlr_u(:) ) * system%Hvol
       enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo  !ilma
    enddo  !ik
   !call comm_summation(u_rVnl_Vnlr_u_l,u_rVnl_Vnlr_u,3*NB*NB*(ik_e-ik_s+1),info%icomm_o)
   !call comm_summation(u_rVnl_Vnlr_u_l,u_rVnl_Vnlr_u,3*NB*NB*NK,info%icomm_ko)


    !(print tm)
    if(flag_print_tm) then

       allocate( u_rVnl_Vnlr_u_all(3,NB,NB,NK), u_rVnl_Vnlr_u_all_l(3,NB,NB,NK) )

       !$omp parallel do private(ik,ib1,ib2) collapse(3)
       do ik =1,NK
       do ib1=1,NB
       do ib2=1,NB
          if( (ik .ge.ik_s .and. ik .le.ik_e) .and.  &
              (ib1.ge.io_s .and. ib1.le.io_e) ) then
             u_rVnl_Vnlr_u_all_l(1:3,ib1,ib2,ik) = u_rVnl_Vnlr_u(1:3,ib1,ib2,ik)
          else
             u_rVnl_Vnlr_u_all_l(1:3,ib1,ib2,ik) = 0d0
          endif
       enddo
       enddo
       enddo
       !$omp end parallel do

       call comm_summation(u_rVnl_Vnlr_u_all_l,u_rVnl_Vnlr_u_all,3*NB*NB*NK,info%icomm_ko)
      !call comm_summation(u_rVnl_Vnlr_u_all_l,u_rVnl_Vnlr_u_all,3*NB*NB*NK,info%icomm_k)

       deallocate(uVpsi_l, uVrpsi_l)
       deallocate(u_rVnl_Vnlr_u_all_l)


       if(yn_ikr=='n')then
          file_tm_data = trim(base_directory)//trim(sysname)//'_tm.data'
       else if(yn_ikr=='y')then
          file_tm_data = trim(base_directory)//trim(sysname)//'_tm_ikr.data'
       end if

       if (comm_is_root(nproc_id_global)) then
          write(*,*) "  printing transition moment ....."

          fh_tm = open_filehandle(file_tm_data, status="replace")
          write(fh_tm, '("#",1X,A)') "#Transition Moment between occupied and unocupied orbitals in GS"
          write(fh_tm, '("#",1X,A)') "# (Separated analysis tool is available)"

          !Currently, TEST: print format is not decided

          !<u_nk|p_j|u_mk>  (j=x,y,z)
          write(fh_tm,*) "#<u_nk|p_j|u_mk>  (j=x,y,z)"
          do ik =1,NK
          do ib1=1,NB
          do ib2=1,NB
             write(fh_tm,9000) ik,ib1,ib2,(upu_all(j,ib1,ib2,ik),j=1,3)
          enddo
          enddo
          enddo
!9000     format(3i8,6e18.10)
9000     format(3i8,6e18.5)

          ! -i*<u_mk|[r_j,V_nl]|u_nk>  (j=x,y,z)
          write(fh_tm,*) "# -i* <u_mk|[r_j,V_nl]|u_nk>  (j=x,y,z)"
          do ik =1,NK
          do ib1=1,NB
          do ib2=1,NB
             write(fh_tm,9000) ik,ib1,ib2,(u_rVnl_Vnlr_u_all(j,ib1,ib2,ik),j=1,3)
          enddo
          enddo
          enddo

          close(fh_tm)
       end if

    end if  !flag_print_tm


    if (flag_print_eps) then
       ! taken from tm2sigma.f90 in utility directory
       if(system%nspin==2) then
          stop "printing option of dielectric function is available for nspin=1"
       endif
       if(xc.ne.'pz' .or. xc.ne.'PZ') then
          if (comm_is_root(nproc_id_global)) then
             write(*,*) "Warning for calculating :"
             write(*,*) "The calculated sigma/epsilon may be wrong except for xc=pz"
             write(*,*) "(not well verified in the current code)"
          endif
       endif


       mu = out_gs_sgm_eps_mu_nu(1)
       nu = out_gs_sgm_eps_mu_nu(2)

       V = system%det_a* dble(NK) ! volume
       n_e = dble(nelec)/system%det_a ! averaged electron number density
  
       if(mu==nu) then
          delta_munu = 1d0
       else
          delta_munu = 0d0
       end if

       !  
       delta  = out_gs_sgm_eps_width
       omega_max = de * nenergy
       nomega = nenergy
       domega = de
       !domega = omega_max/dble(nomega)
       if (comm_is_root(nproc_id_global)) then
       if(mu<1 .or. mu>3 .or. nu<1 .or. nu>3) then
          stop "error: mu & nu must be 1 or 2 or 3 (x,y,z)"
       end if
       end if

       allocate(matrix_vij(io_s:io_e,NB,3,ik_s:ik_e))

       !$omp parallel do private(ik,ib1,ib2) collapse(3)
       do ik=ik_s,ik_e
       do ib1=io_s,io_e
       do ib2=1,NB
          ! <u_ib1,k|p_(1:3)|u_ib2,k> + <u_ib1,k|[r_(1:3),V_nl]|u_ib2,k>/i
          matrix_vij(ib1,ib2,1:3,ik) = upu(1:3,ib1,ib2,ik) + u_rVnl_Vnlr_u(1:3,ib1,ib2,ik)
       end do
       end do
       end do
       !$omp end parallel do
       deallocate( upu, u_rVnl_Vnlr_u )

       allocate(matrix_vji(NB,io_s:io_e,3,ik_s:ik_e))
       allocate(mat_l(NB,3),mat(NB,3))

       do ik=ik_s,ik_e
       do ib2=1,NB

          mat_l(:,:)=0d0
          do ib1=io_s,io_e
             mat_l(ib1,1:3)=matrix_vij(ib1,ib2,1:3,ik)
          enddo

          call comm_summation(mat_l,mat,NB*3,info%icomm_o)

          if(ib2.ge.io_s .and. ib2.le.io_e) then
          do ib1=1,NB
             matrix_vji(ib1,ib2,1:3,ik) = mat(ib1,1:3)
          enddo
          endif

       enddo
       enddo


       if (comm_is_root(nproc_id_global)) then
          filename = trim(base_directory)//trim(sysname) // '_sigma.data'
          open(fh_s, file=filename, status='replace')
          write(fh_s,'(3a)') "#1:omega[a.u.], 2:Re(sigma)[a.u.], 3:Im(sigma)[a.u.]", &
                            & ", 4:Re(sigma_intra)[a.u.], 5:Im(sigma_intra)[a.u.]", &
                            & ", 6:Re(sigma_inter)[a.u.], 7:Im(sigma_inter)[a.u.]"
                            
          filename = trim(base_directory)//trim(sysname) // '_epsilon.data'
          open(fh_e, file=filename, status='replace')
          write(fh_e,'(3a)') "#1:omega[a.u.], 2:Re(epsilon), 3:Im(epsilon)", &
                            & ", 4:Re(epsilon_intra), 5:Im(epsilon_intra)", &
                            & ", 6:Re(epsilon_inter), 7:Im(epsilon_inter)"
       endif

       do iw=1,nomega
          w = dble(iw)*domega
          sigma_l = 0d0
          sigma_intra_l = 0d0
          !$omp parallel do private(ik,ib1,ib2,deigen) collapse(2) reduction(+:sigma_l,sigma_intra_l)
          do ik=ik_s,ik_e
         !do ib1=1,nb
          do ib1=io_s,io_e
             if(system%rocc(ib1,ik,ispin)==0d0) cycle
             do ib2=1,nb
                if(ib2==ib1) cycle
                deigen = energy%esp(ib1,ik,ispin) - energy%esp(ib2,ik,ispin)
                sigma_l = sigma_l + (zi/(w*V))* system%rocc(ib1,ik,ispin)*   &
                     & ( matrix_vij(ib1,ib2,mu,ik) * matrix_vji(ib2,ib1,nu,ik) / ( w + deigen + zi*delta ) &
                     & + matrix_vji(ib2,ib1,mu,ik) * matrix_vij(ib1,ib2,nu,ik) / (-w + deigen - zi*delta ) )
                sigma_intra_l = sigma_intra_l + (zi/(w*V))* system%rocc(ib1,ik,ispin)*   &
                     & ( matrix_vij(ib1,ib2,mu,ik) * matrix_vji(ib2,ib1,nu,ik) / deigen &
                     & + matrix_vji(ib2,ib1,mu,ik) * matrix_vij(ib1,ib2,nu,ik) / deigen )
             end do
          end do
          end do
          !$omp end parallel do
          call comm_summation(sigma_l,      sigma,       info%icomm_ko)
          call comm_summation(sigma_intra_l,sigma_intra, info%icomm_ko)
         !call comm_summation(sigma_l,      sigma,       info%icomm_k)
         !call comm_summation(sigma_intra_l,sigma_intra, info%icomm_k)

          sigma       = sigma       + (zi*n_e/w)* delta_munu
          sigma_intra = sigma_intra + (zi*n_e/w)* delta_munu
          sigma_inter = sigma - sigma_intra
    
          eps = cmplx(delta_munu) + (4d0*pi*zi/w)* sigma
          eps_intra = cmplx(delta_munu) + (4d0*pi*zi/w)* sigma_intra
          eps_inter = (4d0*pi*zi/w)* sigma_inter

          if (comm_is_root(nproc_id_global)) then    
             write(fh_s,'(7(1X,E23.15E3))') w*t_unit_energy%conv,dble(sigma),aimag(sigma), &
                  & dble(sigma_intra),aimag(sigma_intra),dble(sigma_inter),aimag(sigma_inter)
             write(fh_e,'(7(1X,E23.15E3))') w*t_unit_energy%conv,dble(eps),aimag(eps), &
                  & dble(eps_intra),aimag(eps_intra),dble(eps_inter),aimag(eps_inter)
          end if
       end do
  
       close(fh_s)
       close(fh_e)

    endif

    if(yn_ikr=='y')then
       tpsi%zwf(:,:,:,:,:,:,:) = tpsi_backup%zwf(:,:,:,:,:,:,:)
       call deallocate_orbital(tpsi_backup)
    end if

    call comm_sync_all
    return
  end subroutine write_tm_data

  !===================================================================================================================================

  subroutine write_xyz(comment,action,rvf,system,ofl)
  ! Write xyz in xyz format but also velocity and force are printed if necessary
  ! (these can be used for restart of opt and md)
    use structures, only: s_dft_system,s_ofile
    use inputoutput, only: au_length_aa
    use salmon_global, only: base_directory,SYSname,atom_name,base_directory
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem, only: open_filehandle
    implicit none

    type(s_dft_system),intent(in) :: system

    integer :: ia,unit_xyz
    character(3) :: action,rvf
    character(1024) :: file_trj
    character(*) :: comment
    type(s_ofile) :: ofl

    if(.not. comm_is_root(nproc_id_global)) return

    if(action=='new') then
        
       file_trj=trim(base_directory)//trim(SYSname)//'_trj.xyz'
       unit_xyz = open_filehandle(file_trj)
       ofl%fh_trj = unit_xyz
       open(unit_xyz,file=trim(file_trj),status="unknown")

    else if(action=='add') then

       unit_xyz = ofl%fh_trj
       write(unit_xyz,*) system%nion
       write(unit_xyz,*) trim(comment)
       do ia=1,system%nion
          if(      rvf=="r  " ) then
             write(unit_xyz,100) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa
          else if( rvf=="rv " ) then
             write(unit_xyz,110) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa,system%Velocity(1:3,ia)
          else if( rvf=="rvf" ) then
             write(unit_xyz,120) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa,system%Velocity(1:3,ia),system%Force(1:3,ia)
          endif
       enddo

    else if(action=='end') then
       unit_xyz = ofl%fh_trj
       close(unit_xyz)
    endif

100 format(a2,3f18.10)
110 format(a2,3f18.10, "  #v=",3f18.10)
120 format(a2,3f18.10, "  #v=",3f18.10, "  #f=",3f18.10)

  end subroutine write_xyz
  
!===================================================================================================================================

  subroutine write_rt_data_0d(it,ofl,dt,system,rt)
    use structures, only: s_ofile, s_dft_system, s_rt
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem, only: open_filehandle
    use inputoutput, only: t_unit_length,t_unit_time,t_unit_ac,t_unit_elec
    implicit none
    integer, intent(in) :: it
    type(s_ofile) :: ofl
    type(s_dft_system), intent(in) :: system
    type(s_rt),intent(in) :: rt
    integer :: uid, jt
    real(8) :: dm_e(3,system%nspin),ddm_e(3,system%nspin), dm(3), dt

    jt = it
    if(it.lt.0) jt=0
    ddm_e(:,1)= rt%dDp_e(:,jt)
    dm_e(:,1) = rt%Dp_e(:,jt)
    dm(:)     = rt%Dp_e(:,jt) + rt%Dp_i(:,jt) 


    if (comm_is_root(nproc_id_global)) then

    if(it.lt.0) then  !print header
       ofl%fh_rt        = open_filehandle(ofl%file_rt_data)
       uid = ofl%fh_rt

10     format("#",1X,A,":",1X,A)
       write(uid,10) "Real time calculation",""
       write(uid,10) "Ac_ext", "External vector potential field"
       write(uid,10) "E_ext", "External electric field"
       write(uid,10) "Ac_tot", "Total vector potential field"
       write(uid,10) "E_tot", "Total electric field"
       write(uid,10) "ddm_e", "Change of dipole moment (electrons/plus definition)"
       write(uid,10) "dm",    "Total dipole moment (electrons/minus + ions/plus)"

       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 1, "Time", trim(t_unit_time%name), &
         & 2, "Ac_ext_x", trim(t_unit_ac%name), &
         & 3, "Ac_ext_y", trim(t_unit_ac%name), &
         & 4, "Ac_ext_z", trim(t_unit_ac%name), &
         & 5, "E_ext_x", trim(t_unit_elec%name), &
         & 6, "E_ext_y", trim(t_unit_elec%name), &
         & 7, "E_ext_z", trim(t_unit_elec%name), &
         & 8, "Ac_tot_x", trim(t_unit_ac%name), &
         & 9, "Ac_tot_y", trim(t_unit_ac%name), &
         & 10, "Ac_tot_z", trim(t_unit_ac%name), &
         & 11, "E_tot_x", trim(t_unit_elec%name), &
         & 12, "E_tot_y", trim(t_unit_elec%name), &
         & 13, "E_tot_z", trim(t_unit_elec%name), &
         & 14, "ddm_e_x", trim(t_unit_length%name), &
         & 15, "ddm_e_y", trim(t_unit_length%name), &
         & 16, "ddm_e_z", trim(t_unit_length%name), &
         & 17, "dm_x", trim(t_unit_length%name), &
         & 18, "dm_y", trim(t_unit_length%name), &
         & 19, "dm_z", trim(t_unit_length%name)
       write(uid,*)
       flush(uid)

    else  !it>=0
       uid = ofl%fh_rt
       write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & it * dt * t_unit_time%conv,    &
          & system%vec_Ac_ext(1:3) * t_unit_ac%conv, &
          & system%vec_E_ext(1:3) * t_unit_elec%conv, &
          & system%vec_Ac(1:3) * t_unit_ac%conv, &
          & system%vec_E(1:3) * t_unit_elec%conv, &
          & ddm_e(1:3,1) * t_unit_length%conv, &
          & dm(1:3)      * t_unit_length%conv
       write(uid,*)
    endif
    endif

  end subroutine
!===================================================================================================================================
  subroutine write_rt_data_3d(it,ofl,dt,system,curr_e,curr_i)
    use structures, only: s_ofile, s_dft_system
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem, only: open_filehandle
    use inputoutput, only: t_unit_time,t_unit_current,t_unit_ac,t_unit_elec
    use salmon_global, only: spin,yn_md
    implicit none
    type(s_ofile) :: ofl
    integer, intent(in) :: it
    type(s_dft_system), intent(in) :: system
    real(8),intent(in) :: curr_e(3,2), curr_i(3)
    integer :: uid
    real(8) :: dt

    if (comm_is_root(nproc_id_global)) then

    if(it.lt.0) then  !print header
       ofl%fh_rt        = open_filehandle(ofl%file_rt_data)
       uid = ofl%fh_rt

10     format("#",1X,A,":",1X,A)
       write(uid,10) "Real time calculation",""
       write(uid,10) "Ac_ext", "External vector potential field"
       write(uid,10) "E_ext", "External electric field"
       write(uid,10) "Ac_tot", "Total vector potential field"
       write(uid,10) "E_tot", "Total electric field"
       if(spin=='unpolarized' .or. spin=='noncollinear') then
         write(uid,10) "Jm", "Matter current density (electrons)"
       else if(spin=='polarized') then
         write(uid,10) "Jm_u", "Matter current density for spin-up electrons"
         write(uid,10) "Jm_d", "Matter current density for spin-down electrons"
       end if
       if(yn_md=='y') then
         write(uid,10) "Jmi","Matter current density (ions)"
       endif
       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 1, "Time", trim(t_unit_time%name), &
         & 2, "Ac_ext_x", trim(t_unit_ac%name), &
         & 3, "Ac_ext_y", trim(t_unit_ac%name), &
         & 4, "Ac_ext_z", trim(t_unit_ac%name), &
         & 5, "E_ext_x", trim(t_unit_elec%name), &
         & 6, "E_ext_y", trim(t_unit_elec%name), &
         & 7, "E_ext_z", trim(t_unit_elec%name), &
         & 8, "Ac_tot_x", trim(t_unit_ac%name), &
         & 9, "Ac_tot_y", trim(t_unit_ac%name), &
         & 10, "Ac_tot_z", trim(t_unit_ac%name), &
         & 11, "E_tot_x", trim(t_unit_elec%name), &
         & 12, "E_tot_y", trim(t_unit_elec%name), &
         & 13, "E_tot_z", trim(t_unit_elec%name)
       if(spin=='unpolarized' .or. spin=='noncollinear') then
         write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 14, "Jm_x", trim(t_unit_current%name), &
         & 15, "Jm_y", trim(t_unit_current%name), &
         & 16, "Jm_z", trim(t_unit_current%name)
         if(yn_md=='y') then
            write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
                 & 17, "Jmi_x", trim(t_unit_current%name), &
                 & 18, "Jmi_y", trim(t_unit_current%name), &
                 & 19, "Jmi_z", trim(t_unit_current%name)
         endif
       else if(spin=='polarized') then
         write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 14, "Jm_u_x", trim(t_unit_current%name), &
         & 15, "Jm_u_y", trim(t_unit_current%name), &
         & 16, "Jm_u_z", trim(t_unit_current%name), &
         & 17, "Jm_d_x", trim(t_unit_current%name), &
         & 18, "Jm_d_y", trim(t_unit_current%name), &
         & 19, "Jm_d_z", trim(t_unit_current%name)
         if(yn_md=='y') then
            write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
                 & 20, "Jmi_x", trim(t_unit_current%name), &
                 & 21, "Jmi_y", trim(t_unit_current%name), &
                 & 22, "Jmi_z", trim(t_unit_current%name)
         endif
       end if
       write(uid,*)
       flush(uid)

    else  !it>=0
       uid = ofl%fh_rt
       write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & it * dt * t_unit_time%conv,    &
          & system%vec_Ac_ext(1:3) * t_unit_ac%conv, &
          & system%vec_E_ext(1:3) * t_unit_elec%conv, &
          & system%vec_Ac(1:3) * t_unit_ac%conv, &
          & system%vec_E(1:3) * t_unit_elec%conv
       if(spin=='unpolarized') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & curr_e(1:3,1) * t_unit_current%conv
       else if(spin=='polarized') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & curr_e(1:3,1) * t_unit_current%conv, &
          & curr_e(1:3,2) * t_unit_current%conv
       else if(spin=='noncollinear') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & (curr_e(1:3,1)+curr_e(1:3,2)) * t_unit_current%conv
       end if
       if(yn_md=='y') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & curr_i(1:3) * t_unit_current%conv
       endif
       write(uid,*)
    endif
    endif

  end subroutine
  
!===================================================================================================================================

  subroutine write_rt_energy_data(it,ofl,dt,energy,md)
    use structures, only: s_ofile,s_dft_energy,s_md
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use salmon_global, only: ensemble, thermostat, out_rt_energy_step, yn_periodic, yn_jm, yn_fix_func &
    &   ,yn_out_rt_energy_components
    use filesystem, only: open_filehandle
    use inputoutput, only: yn_md,t_unit_time,t_unit_energy
    implicit none
    type(s_dft_energy) :: energy
    type(s_md) :: md
    type(s_ofile) :: ofl
    integer, intent(in) :: it
    integer :: uid,icolumn
    real(8) :: dt

    if(comm_is_root(nproc_id_global)) then

    if(it.lt.0) then  !print header
       ofl%fh_rt_energy = open_filehandle(ofl%file_rt_energy_data)
       uid = ofl%fh_rt_energy

10     format("#",1X,A,":",1X,A)
       write(uid,10) "Real time calculation",""
       write(uid,10) "Eall", "Total energy"
       write(uid,10) "Eall0", "Initial total energy"
       if(yn_out_rt_energy_components=='y') then
         write(uid,10) "E_kin", "Kinetic energy"
         write(uid,10) "E_h", "Hartree energy"
         write(uid,10) "E_ion", "Electron-ion energy"
         write(uid,10) "E_xc", "Exchange-correlation energy"
         write(uid,10) "E_ion_ion", "Ion-ion energy"
       end if
       if(yn_md=='y') then
       write(uid,10) "Tion", "Kinetic energy of ions"
       write(uid,10) "Temperature_ion", "Temperature of ions"
       write(uid,10) "E_work", "Work energy of ions(sum f*dr)"
       if(ensemble=="NVT".and.thermostat=="nose-hoover")then
       write(uid,10) "Enh", "NH thermostat energy (MD)"
       write(uid,10) "Hnvt", "Hamiltonian with NH thermostat(MD)"
       write(uid,10) "Hnvt'","Hnvt using E_work"
       endif
       endif

       if( (yn_periodic=='y' .and. yn_jm=='y') .or. yn_fix_func=='y' ) then
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
          & 1, "Time", trim(t_unit_time%name), &
          & 2, "Eall-Eall0", trim(t_unit_energy%name)
         icolumn=2
       else
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
          & 1, "Time", trim(t_unit_time%name), &
          & 2, "Eall", trim(t_unit_energy%name), &
          & 3, "Eall-Eall0", trim(t_unit_energy%name)
         icolumn=3
       end if

       if(yn_out_rt_energy_components=='y') then
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & icolumn+1, "E_kin", trim(t_unit_energy%name), &
         & icolumn+2, "E_h", trim(t_unit_energy%name), &
         & icolumn+3, "E_ion", trim(t_unit_energy%name), &
         & icolumn+4, "E_xc", trim(t_unit_energy%name), &
         & icolumn+5, "E_ion_ion", trim(t_unit_energy%name)
         icolumn = icolumn + 5
       end if

       if(yn_md=='y') then
       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & icolumn+1, "Tion", trim(t_unit_energy%name), &
        & icolumn+2, "Temperature_ion", "K", &
        & icolumn+3, "E_work", trim(t_unit_energy%name)
       if(ensemble=="NVT".and.thermostat=="nose-hoover")then
       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & icolumn+4, "Enh",  trim(t_unit_energy%name), &
        & icolumn+5, "Hnvt", trim(t_unit_energy%name), &
        & icolumn+6, "Hnvt'",trim(t_unit_energy%name)
       endif
       endif

       write(uid,*)
       flush(uid)
       
       if( (yn_periodic=='y' .and. yn_jm=='y') .or. yn_fix_func=='y' ) then
         write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
             & 0d0,        &
             & 0d0
       else
         write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
             & 0d0,        &
             & energy%E_tot0 * t_unit_energy%conv, &
             & 0d0
       end if
       if(yn_out_rt_energy_components=='y') then
         write(uid, "(99(1X,E23.15E3))",advance='no') &
         & energy%E_kin * t_unit_energy%conv, &
         & energy%E_h * t_unit_energy%conv, &
         & (energy%E_ion_loc+energy%E_ion_nloc) * t_unit_energy%conv, &
         & energy%E_xc * t_unit_energy%conv, &
         & energy%E_ion_ion * t_unit_energy%conv
       end if
       if(yn_md=='y') then
         write(uid, "(99(1X,E23.15E3))",advance='no') &
             & md%Tene * t_unit_energy%conv, &
             & md%Temperature,               &
             & md%E_work * t_unit_energy%conv
       endif
       
       write(uid,*)
       flush(uid)

    else  !it>=0
       if(mod(it,out_rt_energy_step)==0)then
          uid = ofl%fh_rt_energy
   
          if( (yn_periodic=='y' .and. yn_jm=='y') .or. yn_fix_func=='y' ) then
            write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
               & it * dt * t_unit_time%conv,        &
               & (energy%E_tot-energy%E_tot0) * t_unit_energy%conv
          else
            write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
               & it * dt * t_unit_time%conv,        &
               & energy%E_tot * t_unit_energy%conv, &
               & (energy%E_tot-energy%E_tot0) * t_unit_energy%conv
          end if
          if(yn_out_rt_energy_components=='y') then
            write(uid, "(99(1X,E23.15E3))",advance='no') &
            & energy%E_kin * t_unit_energy%conv, &
            & energy%E_h * t_unit_energy%conv, &
            & (energy%E_ion_loc+energy%E_ion_nloc) * t_unit_energy%conv, &
            & energy%E_xc * t_unit_energy%conv, &
            & energy%E_ion_ion * t_unit_energy%conv
          end if
          if(yn_md=='y') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
             & md%Tene * t_unit_energy%conv, &
             & md%Temperature,               &
             & md%E_work * t_unit_energy%conv
          endif
   
   
   !        write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
   !          & it * dt * t_unit_time%conv, &
   !          & Eall_t(it) * t_unit_energy%conv, &
   !          & (Eall_t(it) - Eall0) * t_unit_energy%conv
   !        if(yn_md=='y') then
   !        write(uid, "(99(1X,E23.15E3))",advance='no') &
   !          & Tion_t(it) * t_unit_energy%conv, &
   !          & Temperature_ion_t(it), &
   !          & Ework_integ_fdR(it) * t_unit_energy%conv
   !        if(ensemble=="NVT".and.thermostat=="nose-hoover")then
   !        write(uid, "(99(1X,E23.15E3))",advance='no') &
   !          & Enh_t(it) * t_unit_energy%conv, &
   !          & Hnvt_t(it) * t_unit_energy%conv, &
   !          & (Tion_t(it)+Ework_integ_fdR(it)+Enh_t(it)) * t_unit_energy%conv
   !        endif
   !        endif
          write(uid,*)
          flush(uid)  !for debug
       end if
    endif
    endif

  end subroutine
  
!===================================================================================================================================
  subroutine write_response_0d(ofl,rt)
    use salmon_global, only: e_impulse, nt, dt, nenergy, de
    use inputoutput, only: t_unit_energy,t_unit_polarizability
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz
    real(8) :: tt,hw,t2
    complex(8) :: zalpha(3)
    real(8) :: sf(3)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_response        = open_filehandle(ofl%file_response_data)
      uid = ofl%fh_response

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "alpha", "Polarizability"
      write(uid,10) "df/dE", "Strength function"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Energy", trim(t_unit_energy%name), &
        & 2, "Re(alpha_x)", trim(t_unit_polarizability%name), &
        & 3, "Re(alpha_y)", trim(t_unit_polarizability%name), &
        & 4, "Re(alpha_z)", trim(t_unit_polarizability%name), &
        & 5, "Im(alpha_x)", trim(t_unit_polarizability%name), &
        & 6, "Im(alpha_y)", trim(t_unit_polarizability%name), &
        & 7, "Im(alpha_z)", trim(t_unit_polarizability%name), &
        & 8, "df_x/dE", "none", &
        & 9, "df_y/dE", "none", &
        & 10,"df_z/dE", "none"

      tt = dt*dble(nt)

      do ihw=1,nenergy
        hw=dble(ihw)*de ; zalpha(:)=(0.d0,0.d0)  
        do n=1,nt
          t2=dble(n)*dt ; zalpha(:)=zalpha(:)+exp(zi*hw*t2)*rt%dDp_e(:,n) & 
                                             *(1-3*(t2/tt)**2+2*(t2/tt)**3)
        end do

        zalpha(:)=zalpha(:)/(-e_impulse)*dt
        sf(:)=2*hw/pi*aimag(zalpha(:))

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zalpha(ixyz))*t_unit_polarizability%conv,ixyz=1,3)&
             &,(aimag(zalpha(ixyz))*t_unit_polarizability%conv,ixyz=1,3)&
             &,(sf(ixyz),ixyz=1,3)

      end do

      flush(uid)  !for debug

    end if

  end subroutine

!===================================================================================================================================
  subroutine write_response_3d(ofl,rt)
    use salmon_global, only: e_impulse, trans_longi, nt, dt, nenergy, de, temperature, yn_lr_w0_correction
    use inputoutput, only: t_unit_energy,t_unit_conductivity
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz, it
    real(8) :: tt,hw,t2, smoothing, jav_d(3), jav_s(3), smt_s
    complex(8) :: zsigma(3),zeps(3)
    real(8) :: rtdata(1:3,1:nt)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_response  = open_filehandle(ofl%file_response_data)
      uid = ofl%fh_response

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "sigma", "Conductivity"
      write(uid,10) "eps", "Dielectric constant"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "Energy", trim(t_unit_energy%name), &
        & 2,  "Re(sigma_x)", trim(t_unit_conductivity%name), &
        & 3,  "Re(sigma_y)", trim(t_unit_conductivity%name), &
        & 4,  "Re(sigma_z)", trim(t_unit_conductivity%name), &
        & 5,  "Im(sigma_x)", trim(t_unit_conductivity%name), &
        & 6,  "Im(sigma_y)", trim(t_unit_conductivity%name), &
        & 7,  "Im(sigma_z)", trim(t_unit_conductivity%name), &
        & 8,  "Re(eps_x)", "none", &
        & 9,  "Re(eps_y)", "none", &
        & 10, "Re(eps_z)", "none", &
        & 11, "Im(eps_x)", "none", &
        & 12, "Im(eps_y)", "none", &
        & 13, "Im(eps_z)", "none"

     
      tt = dt*dble(nt)
      if(trans_longi=="tr") then
        rtdata(1:3,1:nt) = rt%curr(1:3,1:nt)
      else if(trans_longi=="lo")then
        rtdata(1:3,1:nt) = rt%E_tot(1:3,1:nt)
      end if

      ! sigma(omega=0) correcton
      jav_d(:)=0d0
      if( yn_lr_w0_correction=='y' .and. temperature < 0d0) then
         jav_s = 0d0
         smt_s = 0d0;
         do it = 1,nt
            t2 = dble(it)*dt
            smoothing = 1d0 - 3d0*(t2/tt)**2 + 2d0*(t2/tt)**3
            jav_s(:)  = jav_s(:) + rtdata(:,it) * smoothing
            smt_s     = smt_s + smoothing
         end do
         jav_d(:) = jav_s(:)/smt_s
      end if

      do ihw=1,nenergy
        hw=dble(ihw)*de
        zsigma(:)=(0.d0,0.d0)
        do n=1,nt
          t2=dble(n)*dt
          smoothing = 1d0 - 3d0*(t2/tt)**2 + 2d0*(t2/tt)**3
         !zsigma(:)=zsigma(:) + exp(zi*hw*t2)* rtdata(:,n) * smoothing
          zsigma(:)=zsigma(:) + exp(zi*hw*t2)* (rtdata(:,n)-jav_d(:)) * smoothing
        end do

        zsigma(:) = (zsigma(:)/e_impulse)*dt
        if(trans_longi=="tr")then
          zeps(:)=1.d0+4.d0*pi*zi*zsigma(:)/hw
        else if(trans_longi=="lo")then
          zeps(:)=1.d0/(1.d0-zsigma(:))
          zsigma(:)=(zeps(:)-1.d0)/(4.d0*pi*zi/hw)
        end if

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zsigma(ixyz))*t_unit_conductivity%conv,ixyz=1,3)&
             &,(aimag(zsigma(ixyz))*t_unit_conductivity%conv,ixyz=1,3)&
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(aimag(zeps(ixyz)),ixyz=1,3)

      end do

      flush(uid)  !for debug

    end if

  end subroutine

!===================================================================================================================================
  subroutine write_dft_md_data(it,ofl,md)
    use structures, only: s_md, s_ofile
    use inputoutput, only: t_unit_time,t_unit_energy
    use salmon_global, only: dt,nt,base_directory,sysname
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root,comm_sync_all
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_md) :: md
    integer :: uid, it

    if(it==0 .and. comm_is_root(nproc_id_global)) then
       ofl%file_dft_md = trim(base_directory)//trim(sysname)//'_dft_md.data'
       ofl%fh_dft_md   = open_filehandle(ofl%file_dft_md)
       uid = ofl%fh_dft_md
       open(uid,file=trim(ofl%file_dft_md),status="unknown")
    endif
    uid = ofl%fh_dft_md

    if(comm_is_root(nproc_id_global)) then
      if(it==0) then
         write(uid,'("#",1X,A)') "DFT-MD: adiabatic (ground state) molecular dynamics"
         write(uid,'("#",1X,A,":",1X,A)') "Tene", "Kinetic energy of atoms(ions)"
         write(uid,'("#",1X,A,":",1X,A)') "Uene", "Potential energy"
         write(uid,'("#",1X,A,":",1X,A)') "Uene0", "Initial potential energy"
         write(uid,'("#",1X,A,":",1X,A)') "E_work", "Work energy (sum F*dR)"
         write(uid,'("#",1X,A,":",1X,A)') "E_nh", "Energy of NH thermostat"
         write(uid,'("#",1X,A,":",1X,A)') "E_tot", "Total energy(=Tene+Uene)"
         write(uid,'("#",1X,A,":",1X,A)') "E_tot0", "Initial total energy"
         write(uid,'("#",1X,A,":",1X,A)') "Hnvt", "Hamiltonian with NH thermostat"
         write(uid,'("#",1X,A,":",1X,A)') "Hnvt'", "Hnvt using E_work"
         write(uid,'("#",1X,A,":",1X,A)') "Temperature", "Temperature of atoms"
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
              &  1, "Time",   trim(t_unit_time%name), &
              &  2, "Tene",   trim(t_unit_energy%name), &
              &  3, "Uene",   trim(t_unit_energy%name), &
              &  4, "Uene-Uene0", trim(t_unit_energy%name), &
              &  5, "E_work", trim(t_unit_energy%name), &
              &  6, "E_nh",   trim(t_unit_energy%name), &
              &  7, "E_tot",  trim(t_unit_energy%name), &
              &  8, "E_tot-E_tot0", trim(t_unit_energy%name), &
              &  9, "Hnvt(2+3+6)",  trim(t_unit_energy%name), &
              & 10, "Hnvt(2+5+6)'", trim(t_unit_energy%name), &
              & 11, "Temperature_ion", "K"
      endif

      write(uid, "(F16.8,99(1X,E23.15E3))") &
          & it * dt    * t_unit_time%conv, &
          & md%Tene    * t_unit_energy%conv, &
          & md%Uene    * t_unit_energy%conv, &
          & (md%Uene-md%Uene0) * t_unit_energy%conv, &
          & md%E_work  * t_unit_energy%conv, &
          & md%E_nh    * t_unit_energy%conv, &
          & md%E_tot   * t_unit_energy%conv, &
          & (md%E_tot-md%E_tot0) * t_unit_energy%conv, &
          & md%Htot    * t_unit_energy%conv, &
          & (md%Tene+md%E_work+md%E_nh) * t_unit_energy%conv, &
          & md%Temperature
      flush(uid)

      if(it==Nt) close(uid)
    end if

    call comm_sync_all
    return
  end subroutine write_dft_md_data
!===================================================================================================================================

  subroutine write_pulse_0d(ofl,rt)
    use inputoutput, only: nt, dt, nenergy, de,  &
                           t_unit_energy,  &
                           t_unit_spectrum_dipole,  &
                           t_unit_spectrum_dipole_square
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz
    real(8) :: tt,hw,t2
    complex(8) :: zdDp_e(3)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_pulse           = open_filehandle(ofl%file_pulse_data)
      uid = ofl%fh_pulse

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "energy", "Frequency"
      write(uid,10) "dm", "Dopile moment"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "energy", trim(t_unit_energy%name), &
        & 2,  "Re(dm_x)", trim(t_unit_spectrum_dipole%name), &
        & 3,  "Re(dm_y)", trim(t_unit_spectrum_dipole%name), &
        & 4,  "Re(dm_z)", trim(t_unit_spectrum_dipole%name), &
        & 5,  "Im(dm_x)", trim(t_unit_spectrum_dipole%name), &
        & 6,  "Im(dm_y)", trim(t_unit_spectrum_dipole%name), &
        & 7,  "Im(dm_z)", trim(t_unit_spectrum_dipole%name), &
        & 8,  "|dm_x|^2", trim(t_unit_spectrum_dipole_square%name), &
        & 9,  "|dm_y|^2", trim(t_unit_spectrum_dipole_square%name), &
        & 10, "|dm_z|^2", trim(t_unit_spectrum_dipole_square%name)

      tt = dt*dble(nt)

      do ihw=1,nenergy
        hw=dble(ihw)*de 
        zdDp_e(:)=(0.d0,0.d0) 
        do n=1,nt
          t2=dble(n)*dt ; zdDp_e(:)=zdDp_e(:)+exp(zi*hw*t2)*rt%dDp_e(:,n) & 
                                             *(1-3*(t2/tt)**2+2*(t2/tt)**3)
        end do
        zdDp_e(:)=zdDp_e(:)*dt

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zdDp_e(ixyz))*t_unit_spectrum_dipole%conv,ixyz=1,3)&
             &,(aimag(zdDp_e(ixyz))*t_unit_spectrum_dipole%conv,ixyz=1,3)&
             &,(abs(zdDp_e(ixyz))**2*t_unit_spectrum_dipole_square%conv,ixyz=1,3)
      end do

    end if

  end subroutine

!===================================================================================================================================
  subroutine write_pulse_3d(ofl,rt)
    use inputoutput, only: nt, dt, nenergy, de,  &
                           t_unit_energy,  &
                           t_unit_spectrum_current,  &
                           t_unit_spectrum_current_square,  &
                           t_unit_spectrum_elec,  &
                           t_unit_spectrum_elec_square
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz
    real(8) :: tt,hw,t2
    complex(8) :: zcurr(3),zE_ext(3),zE_tot(3)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_pulse           = open_filehandle(ofl%file_pulse_data)
      uid = ofl%fh_pulse

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "energy", "Frequency"
      write(uid,10) "Jm", "Matter current"
      write(uid,10) "E_ext", "External electric field"
      write(uid,10) "E_tot", "Total electric field"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "energy", trim(t_unit_energy%name), &
        & 2,  "Re(Jm_x)", trim(t_unit_spectrum_current%name), &
        & 3,  "Re(Jm_y)", trim(t_unit_spectrum_current%name), &
        & 4,  "Re(Jm_z)", trim(t_unit_spectrum_current%name), &
        & 5,  "Im(Jm_x)", trim(t_unit_spectrum_current%name), &
        & 6,  "Im(Jm_y)", trim(t_unit_spectrum_current%name), &
        & 7,  "Im(Jm_z)", trim(t_unit_spectrum_current%name), &
        & 8,  "|Jm_x|^2", trim(t_unit_spectrum_current_square%name), &
        & 9,  "|Jm_y|^2", trim(t_unit_spectrum_current_square%name), &
        & 10, "|Jm_z|^2", trim(t_unit_spectrum_current_square%name), &
        & 11, "Re(E_ext_x)", trim(t_unit_spectrum_elec%name), &
        & 12, "Re(E_ext_y)", trim(t_unit_spectrum_elec%name), &
        & 13, "Re(E_ext_z)", trim(t_unit_spectrum_elec%name), &
        & 14, "Im(E_ext_x)", trim(t_unit_spectrum_elec%name), &
        & 15, "Im(E_ext_y)", trim(t_unit_spectrum_elec%name), &
        & 16, "Im(E_ext_z)", trim(t_unit_spectrum_elec%name), &
        & 17, "|E_ext_x|^2", trim(t_unit_spectrum_elec_square%name), &
        & 18, "|E_ext_y|^2", trim(t_unit_spectrum_elec_square%name), &
        & 19, "|E_ext_z|^2", trim(t_unit_spectrum_elec_square%name), &
        & 20, "Re(E_tot_x)", trim(t_unit_spectrum_elec%name), &
        & 21, "Re(E_tot_y)", trim(t_unit_spectrum_elec%name), &
        & 22, "Re(E_tot_z)", trim(t_unit_spectrum_elec%name), &
        & 23, "Im(E_tot_x)", trim(t_unit_spectrum_elec%name), &
        & 24, "Im(E_tot_y)", trim(t_unit_spectrum_elec%name), &
        & 25, "Im(E_tot_z)", trim(t_unit_spectrum_elec%name), &
        & 26, "|E_tot_x|^2", trim(t_unit_spectrum_elec_square%name), &
        & 27, "|E_tot_y|^2", trim(t_unit_spectrum_elec_square%name), &
        & 28, "|E_tot_z|^2", trim(t_unit_spectrum_elec_square%name)

      tt = dt*dble(nt)

      do ihw=1,nenergy
        hw=dble(ihw)*de 
        zcurr(:)=(0.d0,0.d0) 
        zE_ext(:)=(0.d0,0.d0)  
        zE_tot(:)=(0.d0,0.d0)  
        do n=1,nt
          t2=dble(n)*dt 
          zcurr(:)=zcurr(:)+exp(zi*hw*t2)*rt%curr(:,n) & 
                                      *(1-3*(t2/tt)**2+2*(t2/tt)**3)
          zE_ext(:)=zE_ext(:)+exp(zi*hw*t2)*rt%E_ext(:,n) & 
                                      *(1-3*(t2/tt)**2+2*(t2/tt)**3)
          zE_tot(:)=zE_tot(:)+exp(zi*hw*t2)*rt%E_tot(:,n) & 
                                      *(1-3*(t2/tt)**2+2*(t2/tt)**3)
        end do
        zcurr(:)=zcurr(:)*dt
        zE_ext(:)=zE_ext(:)*dt
        zE_tot(:)=zE_tot(:)*dt

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zcurr(ixyz))*t_unit_spectrum_current%conv,ixyz=1,3)&
             &,(aimag(zcurr(ixyz))*t_unit_spectrum_current%conv,ixyz=1,3)&
             &,(abs(zcurr(ixyz))**2*t_unit_spectrum_current_square%conv,ixyz=1,3)&
             &,(real(zE_ext(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(aimag(zE_ext(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(abs(zE_ext(ixyz))**2*t_unit_spectrum_elec_square%conv,ixyz=1,3)&
             &,(real(zE_tot(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(aimag(zE_tot(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(abs(zE_tot(ixyz))**2*t_unit_spectrum_elec_square%conv,ixyz=1,3)

      end do

      flush(uid)  !for debug

    end if

  end subroutine

!===================================================================================================================================

  subroutine write_prod_dk_data(rgrid_lg, rgrid_mg, system, wf_info, wavefunction)
    use structures,           only: s_rgrid, s_dft_system, s_parallel_info, s_orbital
    use parallelization,      only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem,          only: open_filehandle
    use inputoutput,          only: sysname, base_directory, num_kgrid
    use band,                 only: calc_kgrid_prod
    implicit none
    type(s_rgrid),        intent(in) :: rgrid_lg, rgrid_mg
    type(s_dft_system),       intent(in) :: system
    type(s_parallel_info),      intent(in) :: wf_info
    type(s_orbital), intent(in) :: wavefunction

    ! Specify the neighboring k-grid region to consider:
    integer, parameter :: ndk = 1 
    ! (ndk=1 corresponds to first nearlest neighbors)

    integer :: ik, ik1, ik2, ik3 
    integer :: jdk1, jdk2, jdk3, io, jo
    integer :: fh
    character(256) :: file_prod_dk_data
    integer :: ik3d_tbl(1:3, 1:system%nk)
    complex(8) :: prod_dk( &
      & 1:system%no, 1:system%no, -ndk:ndk, -ndk:ndk, -ndk:ndk, &
      & 1:system%nk)

    ! Export filename: project_directory/sysname_kprod_dk.data
    file_prod_dk_data = trim(base_directory) // trim(sysname) // "_prod_dk.data"

    ! If k-point is distributed as uniform rectangular grid:
    if (0 < minval(num_kgrid)) then
      ! Calculate inner-product table: prod_dk
      call calc_kgrid_prod( &
        & system, rgrid_lg, rgrid_mg, wf_info, wavefunction, &
        & num_kgrid(1), num_kgrid(2), num_kgrid(3), ndk, &
        & ik3d_tbl, prod_dk)
      
      if(comm_is_root(nproc_id_global)) then
        fh = open_filehandle(trim(file_prod_dk_data))
        write(fh, '(a)') "# 1:ik 2:ik1 3:ik2 4:ik3 5:jk1-ik1 6:jk2-ik2 7:jk3-ik3 8:io 9:jo 10:re 11:im"
        do ik = 1, system%nk
          ik1 = ik3d_tbl(1, ik)
          ik2 = ik3d_tbl(2, ik)
          ik3 = ik3d_tbl(3, ik)
          do jdk3 = -ndk, ndk
            do jdk2 = -ndk, ndk
              do jdk1 = -ndk, ndk
                do jo = 1, system%no
                  do io = 1, system%no
                    write(fh, '(9(i10),2(e25.16e3))') &
                      & ik, ik1, ik2, ik3, &
                      & jdk1, jdk2, jdk3, io, jo, &
                      & real(prod_dk(io, jo, jdk1, jdk2, jdk3, ik)), &
                      & aimag(prod_dk(io, jo, jdk1, jdk2, jdk3, ik))
                  end do
                end do
              end do
            end do
          end do
        end do
        close(fh)
      end if
    end if
    return
  end subroutine write_prod_dk_data
  
!===================================================================================================================================

  !! export SYSNAME_info.data file (GS info)
  subroutine write_info_data(Miter,system,energy,pp)
    use structures
    use salmon_global,       only: natom,nelem,iZatom,nelec,sysname,nstate,nelec_spin,unit_system, &
                                   yn_jm, yn_periodic, base_directory
    use parallelization,     only: nproc_id_global
    use communication,only: comm_is_root
    use filesystem,         only: open_filehandle
    use inputoutput,         only: au_length_aa,au_energy_ev
    implicit none
    integer           ,intent(in) :: Miter
    type(s_dft_energy),intent(in) :: energy
    type(s_dft_system),intent(in) :: system
    type(s_pp_info)   ,intent(in) :: pp
    !
    integer :: fh,is,p1,p2,p5,iob,jj,ik,ikoa,iatom,ix
    character(100) :: file_gs_info

    file_gs_info = trim(base_directory)//trim(sysname)//"_info.data"
    fh = open_filehandle(trim(file_gs_info))

    if(comm_is_root(nproc_id_global)) then

       write(fh,*) "Total number of iteration = ", Miter
       write(fh,*)
       write(fh,*) "Number of states = ", nstate
       if(sum(nelec_spin(:))==0) then
          write(fh,*) "Number of electrons = ", nelec
       else
          write(fh,*) "Number of electrons = ", (nelec_spin(is),is=1,2)
       end if
       write(fh,*)
       if(yn_jm=='y'.and.yn_periodic=='y') then
         write(fh,*) "For yn_jm = y and yn_periodic=y, this version still cannot output Total Energy."
       else
         write(fh,*) "Total energy (eV) = ", energy%E_tot*au_energy_ev
       end if
       write(fh,*) "1-particle energies (eV)"
       select case (system%nspin)
       case(1)
          do p5=1,(nstate+3)/4
             p1=4*(p5-1)+1
             p2=4*p5 ; if ( p2 > nstate ) p2=nstate
             write(fh,100) (iob,energy%esp(iob,1,1)*au_energy_ev,iob=p1,p2)
          end do
       case(2)
          write(fh,*) "for up-spin"
          do p5=1,(nstate+3)/4
            p1=4*(p5-1)+1
            p2=4*p5 ; if ( p2 > nstate ) p2=nstate
            write(fh,100) (iob,energy%esp(iob,1,1)*au_energy_ev,iob=p1,p2)
          end do
          write(fh,*) "for down-spin"
          do p5=1,(nstate+3)/4
            p1=4*(p5-1)+1
            p2=4*p5 ; if ( p2 > nstate ) p2=nstate
            write(fh,100) (iob,energy%esp(iob,1,2)*au_energy_ev,iob=p1,p2)
          end do
       end select
       write(fh,*)       
100    format(1x,4(i5,f15.4,2x))

       write(fh,200) "Size of the box (A) = ", system%primitive_a*au_length_aa
       write(fh,200) "Grid spacing (A)    = ", (system%Hgs(jj)*au_length_aa,jj=1,3)
       write(fh,*)
200    format(1x,a,30f14.8)

       if(yn_jm=='n')then
         write(fh,'(1x,"Number of atoms = ",i8)') natom
         do ik=1,nelem
            write(fh,'(1x,"iZatom(",i3,")     = ",i8)') ik, iZatom(ik)
         end do
         write(fh,*)
         write(fh,*) "Ref. and max angular momentum",  &
                     " and pseudo-core radius of PP (A)"
         do ikoa=1,nelem
            write(fh,'(1x,"(",i3,")  "," Ref, Max, Rps =",2i4,f8.3)') &
                  ikoa,pp%Lref(ikoa),pp%Mlps(ikoa),pp%Rps(ikoa)*au_length_aa
         end do
         
         write(fh,*)
         select case(unit_system)
         case('au','a.u.')
            write(fh,*) "Force [au] "
            do iatom=1,natom
               write(fh,300) iatom,(system%Force(ix,iatom),ix=1,3)
            end do
         case('A_eV_fs')
            write(fh,*) "Force [eV/A] "
            do iatom=1,natom
               write(fh,300) iatom,(system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
            end do
         end select
300    format(i6,3e16.8)
       end if

       close(fh)
    endif

  end subroutine write_info_data
  
!===================================================================================================================================

  subroutine write_eigen(ofl,system,energy)
    use structures, only: s_ofile, s_dft_system, s_dft_energy
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use inputoutput, only: uenergy_from_au,iperiodic,unit_energy,sysname,base_directory
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile),intent(inout) :: ofl
    type(s_dft_system),intent(in) :: system
    type(s_dft_energy),intent(in) :: energy
    integer :: iob,iik,is, uid

    if(comm_is_root(nproc_id_global))then
       ofl%file_eigen_data=trim(base_directory)//trim(sysname)//"_eigen.data"
       ofl%fh_eigen = open_filehandle(trim(ofl%file_eigen_data))
       uid = ofl%fh_eigen
       open(uid,file=ofl%file_eigen_data)
       write(uid,'("#esp: single-particle energies (eigen energies)")')
       write(uid,'("#occ: occupation numbers, io: orbital index")')
       select case(unit_energy)
       case('au','a.u.')
          write(uid,'("# 1:io, 2:esp[a.u.], 3:occ")')
       case('ev','eV')
          write(uid,'("# 1:io, 2:esp[eV], 3:occ")')
       end select
       do is=1,system%nspin
       do iik=1,system%nk
          if(iperiodic==3)then
             write(uid,'("k=",1x,i5,",  spin=",1x,i5)') iik,is
          end if
          do iob=1,system%no
             write(uid,'(1x,i5,2(e26.16e3))') iob, energy%esp(iob,iik,is)*uenergy_from_au, system%rocc(iob,iik,is)
          end do
       end do
       end do
       close(uid)
    end if

  end subroutine write_eigen
  
!===================================================================================================================================

  subroutine write_dos(system,energy)
    use structures
    use math_constants, only: pi
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use inputoutput, only: uenergy_from_au
    use salmon_global, only: out_dos_start, out_dos_end, out_dos_function, &
                           out_dos_width, out_dos_nenergy, yn_out_dos_set_fe_origin, unit_energy, &
                           nelec,nstate,temperature,yn_spinorbit, base_directory,sysname
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_dft_energy),intent(in) :: energy
    !
    integer :: iob,iik,is
    real(8) :: dos_l(1:out_dos_nenergy,system%nspin)
    real(8) :: fk,ww,dw
    integer :: iw,index_vbm
    real(8) :: ene_min,ene_max,eshift
    character(100) :: filename

    ene_min = minval(energy%esp)
    ene_max = maxval(energy%esp)
    if(yn_out_dos_set_fe_origin=='y') then
      if(temperature>=0.d0) then
        eshift = system%mu
      else
        if( yn_spinorbit=='y' )then
          index_vbm=nelec
        else
          index_vbm=nelec/2
        end if
        eshift = maxval(energy%esp(1:index_vbm,:,:)) ! valence band maximum (HOMO)
      end if
    else
      eshift = 0.d0
    endif
    out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min)-eshift)
    out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min)-eshift)
    dw=(out_dos_end-out_dos_start)/dble(out_dos_nenergy-1)

    dos_l = 0.d0

    do is=1,system%nspin
    do iik=1,system%nk
    do iob=1,system%no
      select case (out_dos_function)
      case('lorentzian')
        fk=2.d0*out_dos_width/pi
        do iw=1,out_dos_nenergy
          ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,is)
          dos_l(iw,is)=dos_l(iw,is)+system%wtk(iik)*fk/(ww**2+out_dos_width**2)
        end do
      case('gaussian')
        fk=2.d0/(sqrt(2.d0*pi)*out_dos_width)
        do iw=1,out_dos_nenergy
          ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,is)
          dos_l(iw,is)=dos_l(iw,is)+system%wtk(iik)*fk*exp(-(0.5d0/out_dos_width**2)*ww**2)
        end do
      end select
    end do
    end do
    end do

    if(comm_is_root(nproc_id_global))then
      filename=trim(base_directory)//trim(sysname)//"_dos.data"
      open(101,file=filename)
      write(101,'("# Density of States")')
      select case(unit_energy)
      case('au','a.u.')
        write(101,'("# Energy[a.u.] DOS[a.u.]")')
      case('ev','eV')
        write(101,'("# Energy[eV]  DOS[1/eV]")')
      end select
      write(101,'("#-----------------------")')
      do iw=1,out_dos_nenergy
        ww=out_dos_start+dble(iw-1)*dw
        write(101,'(F16.8,99(1X,E23.15E3))') ww*uenergy_from_au, ( dos_l(iw,is)/uenergy_from_au, is=1,system%nspin )
      end do
      close(101)
    end if

  end subroutine write_dos
  
!===================================================================================================================================
  
  subroutine write_pdos(lg,mg,system,info,pp,energy,tpsi)
    use structures
    use math_constants      ,only: pi,zi
    use salmon_math         ,only: ylm
    use parallelization     ,only: nproc_id_global
    use communication       ,only: comm_is_root, comm_summation
    use salmon_global       ,only: out_dos_start, out_dos_end, out_dos_function, &
                                   out_dos_width, out_dos_nenergy, yn_out_dos_set_fe_origin, &
                                   nelec, kion, natom, nstate, unit_energy, temperature, yn_spinorbit, &
                                   base_directory,sysname
    use inputoutput         ,only: uenergy_from_au
    use prep_pp_sub         ,only: bisection
    implicit none
    type(s_rgrid)           ,intent(in) :: lg,mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_pp_info)         ,intent(in) :: pp
    type(s_dft_energy)      ,intent(in) :: energy
    type(s_orbital)         ,intent(in) :: tpsi
    !
    integer :: iob,iatom,L,m,ix,iy,iz,iik,ispin
    integer :: ikoa
    integer :: intr
    real(8) :: phi_r
    real(8) :: rr
    real(8) :: ratio1,ratio2
    real(8) :: xx,yy,zz
    real(8) :: xxxx,yyyy,zzzz,rinv
    integer :: lm
    real(8) :: rbox_pdos(25,natom)
    real(8) :: rbox_pdos2(25,natom)
    real(8) :: pdos_l_tmp(out_dos_nenergy,0:4,natom)
    real(8) :: pdos_l(out_dos_nenergy,0:4,natom)
    character(100) :: Outfile
    real(8) :: fk,ww,dw
    integer :: iw,index_vbm
    real(8) :: ene_min,ene_max,eshift
    character(20) :: fileNumber

    if( all(pp%upp_f==0.0d0) )then
      write(*,*) "@calc_pdos: Pseudoatom wave function is not available"
      return
    end if

    ene_min = minval(energy%esp(:,:,:))
    ene_max = maxval(energy%esp(:,:,:))
    if(yn_out_dos_set_fe_origin=='y') then
      if(temperature>=0.d0) then
        eshift = system%mu
      else
        if( yn_spinorbit=='y' )then
          index_vbm=nelec
        else
          index_vbm=nelec/2
        end if
        eshift = maxval(energy%esp(1:index_vbm,:,:)) ! valence band maximum (HOMO)
      end if
    else
      eshift = 0.d0
    endif
    out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min)-eshift)
    out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min)-eshift)
    dw=(out_dos_end-out_dos_start)/dble(out_dos_nenergy-1)

    pdos_l_tmp=0.d0

    do ispin=1,system%nspin
    do iik=info%ik_s,info%ik_e
    do iob=info%io_s,info%io_e
      rbox_pdos=0.d0
      do iatom=1,natom
        ikoa=Kion(iatom)
        do L=0,pp%mlps(ikoa)
          do m=-L,L
            lm=L*L+L+1+m
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              xx=lg%coordinate(ix,1)-system%Rion(1,iatom)
              yy=lg%coordinate(iy,2)-system%Rion(2,iatom)
              zz=lg%coordinate(iz,3)-system%Rion(3,iatom)
              rr=sqrt(xx**2+yy**2+zz**2)+1.d-50
              rinv=1.0d0/rr
              xxxx=xx*rinv
              yyyy=yy*rinv
              zzzz=zz*rinv
              call bisection(rr,intr,ikoa,pp%nrmax,pp%rad)
              if(intr==1) intr=2
              ratio1=(rr-pp%rad(intr,ikoa))/(pp%rad(intr+1,ikoa)-pp%rad(intr,ikoa)) ; ratio2=1.d0-ratio1
              phi_r= ratio1*pp%upp_f(intr,pp%lref(ikoa),ikoa)/rr**(pp%lref(ikoa)+1)*sqrt((2*pp%lref(ikoa)+1)/(4*Pi)) +  &
                     ratio2*pp%upp_f(intr-1,pp%lref(ikoa),ikoa)/rr**(pp%lref(ikoa)+1)*sqrt((2*pp%lref(ikoa)+1)/(4*Pi))
                                            !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
              if(allocated(tpsi%rwf)) then
                rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+tpsi%rwf(ix,iy,iz,ispin,iob,iik,1)*phi_r*Ylm(xxxx,yyyy,zzzz,L,m)*system%Hvol
              else
                rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+tpsi%zwf(ix,iy,iz,ispin,iob,iik,1)*phi_r*Ylm(xxxx,yyyy,zzzz,L,m)*system%Hvol
              end if
            end do
            end do
            end do
          end do
        end do
      end do
      call comm_summation(rbox_pdos,rbox_pdos2,25*natom,info%icomm_r)
      do iatom=1,natom
        ikoa=Kion(iatom)
        do L=0,pp%mlps(ikoa)
          do lm=L**2+1,(L+1)**2
            select case (out_dos_function)
            case('lorentzian')
              fk=2.d0*out_dos_width/pi
              do iw=1,out_dos_nenergy
                ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,ispin)
                pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
                  +abs(rbox_pdos2(lm,iatom))**2*fk/(ww**2+out_dos_width**2)
              end do
            case('gaussian')
              fk=2.d0/(sqrt(2.d0*pi)*out_dos_width)
              do iw=1,out_dos_nenergy
                ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,ispin)
                pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
                  +abs(rbox_pdos2(lm,iatom))**2*fk*exp(-(0.5d0/out_dos_width**2)*ww**2)
              end do
            end select
          end do
        end do
      end do
    end do
    end do
    end do
    call comm_summation(pdos_l_tmp,pdos_l,out_dos_nenergy*5*natom,info%icomm_ko)

    if(comm_is_root(nproc_id_global))then
      do iatom=1,natom
        ikoa=Kion(iatom)
        write(fileNumber, '(i8)') iatom
        OutFile = trim(base_directory)//trim(sysname)//"_pdos"//trim(adjustl(fileNumber))//".data"
        open(101,file=OutFile)
        write(101,'("# Projected Density of States")')
        select case(unit_energy)
        case('au','a.u.')
          if(pp%mlps(ikoa)==0)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.]")')
          else if(pp%mlps(ikoa)==1)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.]")')
          else if(pp%mlps(ikoa)==2)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.]")')
          else if(pp%mlps(ikoa)==3)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.] PDOS(l=3)[a.u.]")')
          end if
        case('ev','eV')
          if(pp%mlps(ikoa)==0)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV]")')
          else if(pp%mlps(ikoa)==1)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV]")')
          else if(pp%mlps(ikoa)==2)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV]")')
          else if(pp%mlps(ikoa)==3)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV] PDOS(l=3)[1/eV]")')
          end if
        end select
        write(101,'("#-----------------------")')
        if(pp%mlps(ikoa)==0)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        else if(pp%mlps(ikoa)==1)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,2f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        else if(pp%mlps(ikoa)==2)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,3f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        else if(pp%mlps(ikoa)==3)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,4f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        end if
        close(101)
      end do
    end if
    return
  end subroutine write_pdos
  
!===================================================================================================================================

  subroutine write_band_information(system,energy)
    use structures
    use salmon_global, only: nelec,yn_periodic,yn_spinorbit
    use inputoutput, only: au_energy_ev
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_dft_energy),intent(in) :: energy
    !
    integer :: ik,index_vbm
    real(8),dimension(system%nk) :: esp_vb_min,esp_vb_max,esp_cb_min,esp_cb_max
    !
    if( yn_spinorbit=='y' )then
      index_vbm=nelec
    else
      index_vbm=nelec/2
      if(mod(nelec,2)==1) index_vbm = index_vbm + 1
    end if
    if(comm_is_root(nproc_id_global) .and. index_vbm<system%no) then
      do ik=1,system%nk
        esp_vb_min(ik)=minval(energy%esp(1:index_vbm,ik,:))
        esp_vb_max(ik)=maxval(energy%esp(1:index_vbm,ik,:))
        esp_cb_min(ik)=minval(energy%esp(index_vbm+1:system%no,ik,:))
        esp_cb_max(ik)=maxval(energy%esp(index_vbm+1:system%no,ik,:))
      end do
      if(yn_periodic=='y') then
        write(*,*) 'band information-----------------------------------------'
        write(*,*) 'Bottom of VB',minval(esp_vb_min(:))
        write(*,*) 'Top of VB',maxval(esp_vb_max(:))
        write(*,*) 'Bottom of CB',minval(esp_cb_min(:))
        write(*,*) 'Top of CB',maxval(esp_cb_max(:))
        write(*,*) 'Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
        write(*,*) 'BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
        write(*,*) 'Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
        write(*,*) 'Physicaly upper bound of eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
        write(*,*) '---------------------------------------------------------'
        write(*,*) 'Bottom of VB[eV]',minval(esp_vb_min(:))*au_energy_ev
        write(*,*) 'Top of VB[eV]',maxval(esp_vb_max(:))*au_energy_ev
        write(*,*) 'Bottom of CB[eV]',minval(esp_cb_min(:))*au_energy_ev
        write(*,*) 'Top of CB[eV]',maxval(esp_cb_max(:))*au_energy_ev
        write(*,*) 'Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*au_energy_ev
        write(*,*) 'BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*au_energy_ev
        write(*,*) 'Physicaly upper bound of CB for DOS[eV]',minval(esp_cb_max(:))*au_energy_ev
        write(*,*) 'Physicaly upper bound of eps(omega)[eV]',minval(esp_cb_max(:)-esp_vb_min(:))*au_energy_ev
        write(*,*) '---------------------------------------------------------'
      else
        if(system%nk /= 1) stop "error: yn_periodic='n' and Nk/=1"
        write(*,*) 'orbital energy information-------------------------------'
        write(*,*) 'Lowest occupied orbital',esp_vb_min(1)
        write(*,*) 'Highest occupied orbital (HOMO)',esp_vb_max(1)
        write(*,*) 'Lowest unoccupied orbital (LUMO)',esp_cb_min(1)
        write(*,*) 'Highest unoccupied orbital',esp_cb_max(1)
        write(*,*) 'HOMO-LUMO gap',esp_cb_min(1)-esp_vb_max(1)
        write(*,*) 'Physicaly upper bound of eps(omega)',esp_cb_max(1)-esp_vb_min(1)
        write(*,*) '---------------------------------------------------------'
        write(*,*) 'Lowest occupied orbital[eV]',esp_vb_min(1)*au_energy_ev
        write(*,*) 'Highest occupied orbital (HOMO)[eV]',esp_vb_max(1)*au_energy_ev
        write(*,*) 'Lowest unoccupied orbital (LUMO)[eV]',esp_cb_min(1)*au_energy_ev
        write(*,*) 'Highest unoccupied orbital[eV]',esp_cb_max(1)*au_energy_ev
        write(*,*) 'HOMO-LUMO gap[eV]',(esp_cb_min(1)-esp_vb_max(1))*au_energy_ev
        write(*,*) 'Physicaly upper bound of eps(omega)[eV]',(esp_cb_max(1)-esp_vb_min(1))*au_energy_ev
        write(*,*) '---------------------------------------------------------'
      end if
    end if
    return
  end subroutine write_band_information
  
!===================================================================================================================================
  subroutine init_projection(system,lg,mg,info,stencil,Vpsl,xc_func,pp,ppn,fg,poisson,srg_scalar,rt,energy,ofl)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: projection_option,nstate,directory_read_data,yn_restart,yn_self_checkpoint
    use checkpoint_restart_sub, only: read_bin,generate_restart_directory_name
    use initialization_sub, only: init_parallel_dft
    use sendrecv_grid, only: create_sendrecv_neig, init_sendrecv_grid
    implicit none
    type(s_rgrid)           ,intent(in) :: lg,mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_scalar)          ,intent(in) :: Vpsl
    type(s_xc_functional)   ,intent(in) :: xc_func
    type(s_pp_info)         ,intent(in) :: pp
    type(s_pp_nlcc)         ,intent(in) :: ppn
    type(s_reciprocal_grid) ,intent(in) :: fg
    type(s_poisson)                     :: poisson
    type(s_sendrecv_grid)               :: srg_scalar
    type(s_rt)                          :: rt
    type(s_dft_energy)                  :: energy
    type(s_ofile)                       :: ofl
    !
    character(256) :: wdir,gdir,dir_gs
    logical :: iself
    integer :: iter,jspin
    integer,dimension(2,3) :: neig
    
    call file_header
    
    if(yn_restart=='n') then
      dir_gs = directory_read_data
    else
      dir_gs = 'restart/' ! default GS directory
      if(comm_is_root(nproc_id_global)) write(*,*) " projection_option: read GS data from directory ./restart/"
    end if
    
    rt%system_proj = system
    rt%system_proj%no = nstate ! # of orbitals for GS ! future work: nstate --> nstate_proj
    deallocate(rt%system_proj%rocc)
    allocate(rt%system_proj%rocc(rt%system_proj%no,system%nk,system%nspin)) ! this will be filled by read_bin
    
    rt%info_proj = info
    deallocate(rt%info_proj%io_s_all,rt%info_proj%io_e_all,rt%info_proj%numo_all,rt%info_proj%irank_io)
    call init_parallel_dft(rt%system_proj,rt%info_proj)
    
    call create_sendrecv_neig(neig, rt%info_proj)
    call init_sendrecv_grid(rt%srg_proj, mg, &
    & rt%info_proj%numo*rt%info_proj%numk*rt%system_proj%nspin, rt%info_proj%icomm_rko, neig)

    call allocate_orbital_complex(system%nspin,mg,rt%info_proj,rt%tpsi0)
    call allocate_orbital_complex(system%nspin,mg,rt%info_proj,rt%ttpsi0)
    call allocate_orbital_complex(system%nspin,mg,rt%info_proj,rt%htpsi0)

    deallocate(energy%esp)
    allocate(energy%esp(rt%system_proj%no,system%nk,system%nspin))
    
  ! wavefunctions @ GS
    call generate_restart_directory_name(dir_gs,gdir,wdir)
    iself = yn_restart =='y' .and. yn_self_checkpoint == 'y'
    if (.not. iself) then
      wdir = gdir
    end if
    call read_bin(wdir,lg,mg,rt%system_proj,rt%info_proj,rt%tpsi0,iter,is_self_checkpoint=iself)

  ! V_local
    allocate(rt%vloc0(system%nspin))
    do jspin=1,system%nspin
      call allocate_scalar(mg,rt%vloc0(jspin))
    end do
    if(projection_option=='gs') then
      call calc_vloc0 ! V_local for GS
    end if
    
  contains
  
    subroutine file_header
      use salmon_global, only: base_directory,SYSname,yn_out_intraband_current,spin
      use inputoutput, only: t_unit_time,t_unit_current
      use filesystem, only: open_filehandle
      implicit none
      character(256) :: file_intra_current
      if(comm_is_root(nproc_id_global)) then
      
      !(header in SYSname_ovlp.data)
        write(ofl%file_ovlp,"(2A,'_ovlp.data')") trim(base_directory),trim(SYSname)
        ofl%fh_ovlp = open_filehandle(ofl%file_ovlp)
        open(ofl%fh_ovlp,file=ofl%file_ovlp)
        write(ofl%fh_ovlp, '("#",1X,A)') "Projection"
        write(ofl%fh_ovlp, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(ofl%fh_ovlp, '("#",1X,A,":",1X,A)') "ovlp_occup", "Occupation"
        write(ofl%fh_ovlp, '("#",1X,A,":",1X,A)') "NB", "Number of bands"
        write(ofl%fh_ovlp, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "ik", "none", &
        & 2, "ovlp_occup(1:NB)", "none"
        
      !(header in SYSname_nex.data)
        write(ofl%file_nex,"(2A,'_nex.data')") trim(base_directory),trim(SYSname)
        ofl%fh_nex = open_filehandle(ofl%file_nex)
        open(ofl%fh_nex,file=ofl%file_nex)
        write(ofl%fh_nex, '("#",1X,A)') "Excitation"
        write(ofl%fh_nex, '("#",1X,A,":",1X,A)') "nelec", "Number of excited electrons"
        write(ofl%fh_nex, '("#",1X,A,":",1X,A)') "nhole", "Number of excited holes"
        write(ofl%fh_nex, '("#",99(1X,I0,":",A,"[",A,"]"))')  &
        &           1, "time", trim(t_unit_time%name), &
        &           2, "nelec", "none", &
        &           3, "nhole", "none"
        
        if(yn_out_intraband_current=='y') then
        !(header in SYSname_intra_current.data)
          write(file_intra_current,"(2A,'_intra_current.data')") trim(base_directory),trim(SYSname)
          ofl%fh_intra_current = open_filehandle(file_intra_current)
          open(ofl%fh_intra_current,file=file_intra_current)
          write(ofl%fh_intra_current, '("#",1X,A)') "Intra-band current density [T. Otobe, Phys. Rev. B 94, 235152 (2016).]"
          write(ofl%fh_intra_current, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no')  &
          &  1, "Time", trim(t_unit_time%name)
          if(spin=='unpolarized' .or. spin=='noncollinear') then
            write(ofl%fh_intra_current, '(99(1X,I0,":",A,"[",A,"]"))') &
            & 2, "Jm_intra_x", trim(t_unit_current%name), &
            & 3, "Jm_intra_y", trim(t_unit_current%name), &
            & 4, "Jm_intra_z", trim(t_unit_current%name)
          else if(spin=='polarized') then
            write(ofl%fh_intra_current, '(99(1X,I0,":",A,"[",A,"]"))') &
            & 2, "Jm_intra_u_x", trim(t_unit_current%name), &
            & 3, "Jm_intra_u_y", trim(t_unit_current%name), &
            & 4, "Jm_intra_u_z", trim(t_unit_current%name), &
            & 5, "Jm_intra_d_x", trim(t_unit_current%name), &
            & 6, "Jm_intra_d_y", trim(t_unit_current%name), &
            & 7, "Jm_intra_d_z", trim(t_unit_current%name)
          end if
        end if
        
      end if
    end subroutine file_header
  
    subroutine calc_vloc0
      use density_matrix, only: calc_density
      use hartree_sub, only: hartree
      use salmon_xc, only: exchange_correlation
      use hamiltonian, only: update_vlocal
      implicit none
      real(8) :: E_xc
      type(s_scalar) :: rho,Vh
      type(s_scalar),allocatable :: rho_s(:),Vxc(:)
      
      call allocate_scalar(mg,rho)
      call allocate_scalar(mg,Vh)
      allocate(rho_s(system%nspin),Vxc(system%nspin))
      do jspin=1,system%nspin
        call allocate_scalar(mg,rho_s(jspin))
        call allocate_scalar(mg,Vxc(jspin))
      end do
      
      call calc_density(rt%system_proj,rho_s,rt%tpsi0,rt%info_proj,mg)
      rho%f = 0d0
      do jspin=1,system%nspin
        rho%f = rho%f + rho_s(jspin)%f
      end do
      call hartree(lg,mg,rt%info_proj,rt%system_proj,fg,poisson,srg_scalar,stencil,rho,Vh)
      call exchange_correlation(rt%system_proj,xc_func,mg,srg_scalar,rt%srg_proj,rho_s &
      & ,pp,ppn,rt%info_proj,rt%tpsi0,stencil,Vxc,E_xc)
      call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,rt%vloc0)
      
      call deallocate_scalar(rho)
      call deallocate_scalar(Vh)
      do jspin=1,system%nspin
        call deallocate_scalar(rho_s(jspin))
        call deallocate_scalar(Vxc(jspin))
      end do
      deallocate(rho_s,Vxc)
    end subroutine calc_vloc0

  end subroutine init_projection

  subroutine projection(itt,ofl,dt,mg,system,info,stencil,V_local,ppg,psi_t,energy,rt)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: projection_option,ncg,nelec,yn_spinorbit,nscf,yn_out_intraband_current,threshold_projection
    use inputoutput, only: t_unit_time
    use subspace_diagonalization, only: ssdg
    use gram_schmidt_orth, only: gram_schmidt
    use Conjugate_Gradient, only: gscg_zwf
    use Total_Energy, only: calc_eigen_energy
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)           ,intent(in) :: ofl
    real(8)                 ,intent(in) :: dt
    type(s_rgrid)           ,intent(in) :: mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_scalar)          ,intent(in) :: V_local(system%nspin)
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_orbital)         ,intent(in) :: psi_t ! | u_{n,k}(t) >
    type(s_dft_energy)                  :: energy
    type(s_rt)                          :: rt
    !
    integer :: nspin,nspin_tmp,no,no0,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
    integer :: ix,iy,iz,io1,io2,io,ik,ispin,iter_GS,niter
    complex(8),dimension(rt%system_proj%no,system%no,system%nspin,system%nk) :: mat
    real(8) :: coef(rt%system_proj%no,system%nk,system%nspin)
    real(8) :: nee, neh, wspin, dE
    complex(8) :: cbox
      
    if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ projection"
    
    rt%system_proj%vec_Ac(1:3) = system%vec_Ac(1:3) ! vector potential @ time itt

    no0 = rt%system_proj%no ! # of orbitals @ GS
    
    nspin = system%nspin
    no = system%no
    nk = system%nk
    is = mg%is
    ie = mg%ie
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    
    if(nspin==1) then
      wspin = 2d0
    else if(nspin==2) then
      wspin = 1d0
    end if
    
    nspin_tmp = nspin
    if(yn_spinorbit=='y') then
      nspin_tmp = 1
    end if
    
    if(nscf==0) then
      niter = 10
    else
      niter = nscf
    end if

    if(projection_option=='td') then
      do ispin=1,nspin
        rt%vloc0(ispin)%f = V_local(ispin)%f ! the local potential for the time-dependent Hamiltonian
      end do
    end if

!
! psi_t    = | u_{n,k}(t) >, TDDFT orbitals
!
! rt%tpsi0 = | w_{n,k}(t) >, projection orbitals
!
! projection_option=='gs':
!   | w_{n,k}(t) > = the eigenorbitals for the ground-state Hamiltonian whose k-point is shifted by A(t).
!   i.e. | w_{n,k}(t) > = | u^{GS}_{n,k+A(t)/c} >
!   a.k.a the Houston functions
!
! projection_option=='td':
!   | w_{n,k}(t) > = the instantaneous eigenorbitals for the time-dependent Hamiltonian.
!


    call calc_eigen_energy(energy,rt%tpsi0,rt%htpsi0,rt%ttpsi0 &
    & ,rt%system_proj,rt%info_proj,mg,rt%vloc0,stencil,rt%srg_proj,ppg)
    dE = energy%E_kin - rt%E_old
    if(abs(dE) < min(threshold_projection,1e-12)) then
      if(comm_is_root(nproc_id_global)) write(*,*) "projection: already converged, E_kin(new)-E_kin(old)=",dE
    else
      do iter_GS=1,niter
        call ssdg(mg,rt%system_proj,rt%info_proj,stencil,rt%tpsi0,rt%htpsi0,ppg,rt%vloc0,rt%srg_proj)
        call gscg_zwf(ncg,mg,rt%system_proj,rt%info_proj,stencil,ppg,rt%vloc0,rt%srg_proj,rt%tpsi0,rt%cg)
        call gram_schmidt(rt%system_proj, mg, rt%info_proj, rt%tpsi0)
        call calc_eigen_energy(energy,rt%tpsi0,rt%htpsi0,rt%ttpsi0 &
        & ,rt%system_proj,rt%info_proj,mg,rt%vloc0,stencil,rt%srg_proj,ppg)
        dE = energy%E_kin - rt%E_old
        if(comm_is_root(nproc_id_global)) write(*,'(a,i6,e20.10)') "projection: ",iter_GS,dE
        if(abs(dE) < threshold_projection) exit
        rt%E_old = energy%E_kin
      end do
    end if
    
    call inner_product(rt%tpsi0,psi_t,mat) ! mat(n,m) = < w_{n,k)(t) | u_{m,k}(t) >

    if(yn_spinorbit=='y') then
      mat(1:no0,1:no,1,1:nk) = mat(1:no0,1:no,1,1:nk) + mat(1:no0,1:no,2,1:nk)
      mat(1:no0,1:no,2,1:nk) = mat(1:no0,1:no,1,1:nk)
    end if

    coef=0.d0
    do ispin=1,nspin
    do ik=1,nk
    do io1=1,no0
      do io2=1,no
        coef(io1,ik,ispin) = coef(io1,ik,ispin) &
        & + system%rocc(io2,ik,ispin)*system%wtk(ik)* abs(mat(io1,io2,ispin,ik))**2
      end do
    end do
    end do
    end do

    nee = 0d0
    neh = dble(nelec)
    do ispin=1,nspin_tmp
    do ik=1,nk
    do io=1,no0
    ! /wspin: for canceling double counting of 2 in rocc.
      nee = nee + ((wspin-rt%system_proj%rocc(io,ik,ispin))/wspin) * coef(io,ik,ispin)
      neh = neh - (rt%system_proj%rocc(io,ik,ispin)/wspin) * coef(io,ik,ispin)
    end do
    end do
    end do
   !nee  = sum(ovlp_occ(NBoccmax+1:NB,:))
   !neh  = sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
    if(comm_is_root(nproc_id_global))then
      write(ofl%fh_nex,'(99(1X,E23.15E3))') dble(itt)*dt*t_unit_time%conv, nee, neh
      write(ofl%fh_ovlp,'(i11)') itt
      do ispin=1,nspin_tmp
      do ik=1,nk
        write(ofl%fh_ovlp,'(i6,1000(1X,E23.15E3))') ik,(coef(io,ik,ispin)/system%wtk(ik),io=1,no0)
      end do
      end do
    end if
    
    if(yn_out_intraband_current=='y') then
      call write_intraband_current
    end if
    
!        if(action=="proj_last ") then
!       tconv = t_unit_energy%conv
!       do ik=1,NK
!          write(409,10)ik,(esp_all(ib,ik)*tconv,ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
!       end do

    return
    
  contains
  
    subroutine inner_product(psi1,psi2,mat)
      use communication, only: comm_summation, comm_bcast
      use pack_unpack, only: copy_data
      implicit none
      type(s_orbital), intent(in) :: psi1,psi2
      complex(8),dimension(rt%system_proj%no,system%no,system%nspin,system%nk) :: mat
      !
      complex(8),dimension(rt%system_proj%no,system%no,system%nspin,system%nk) :: mat1
      complex(8) :: wf_io1(mg%is_array(1):mg%ie_array(1) &
                        & ,mg%is_array(2):mg%ie_array(2) &
                        & ,mg%is_array(3):mg%ie_array(3))
      ! copied from subspace_diagonalization.f90
      mat1 = 0d0
      if(info%if_divide_orbit) then
        do ik=ik_s,ik_e
        do ispin = 1, nspin
          do io1 = 1, no0
            if (rt%info_proj%io_s<= io1 .and. io1 <= rt%info_proj%io_e) then
              call copy_data(psi1%zwf(:, :, :, ispin, io1, ik, 1),wf_io1)
            end if
            call comm_bcast(wf_io1, rt%info_proj%icomm_o, rt%info_proj%irank_io(io1))
            do io2 = 1, no
              if (io_s<= io2 .and. io2 <= io_e) then
                cbox = 0d0
                !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+:cbox)
                do iz=is(3),ie(3)
                do iy=is(2),ie(2)
                do ix=is(1),ie(1)
                  cbox = cbox + conjg(wf_io1(ix,iy,iz)) * psi2%zwf(ix,iy,iz,ispin,io2,ik,1)
                end do
                end do
                end do
                mat1(io1,io2,ispin,ik) = cbox * system%hvol
              end if
            end do
          end do !io1
        end do !ispin
        end do
      else
        !$omp parallel do private(ik,io1,io2,ispin,cbox,iz,iy,ix) collapse(4)
        do ik=ik_s,ik_e
        do ispin=1,nspin
        do io1=1,no0
        do io2=1,no
          cbox = 0d0
          do iz=is(3),ie(3)
          do iy=is(2),ie(2)
          do ix=is(1),ie(1)
            cbox = cbox + conjg(psi1%zwf(ix,iy,iz,ispin,io1,ik,1)) * psi2%zwf(ix,iy,iz,ispin,io2,ik,1)
          end do
          end do
          end do
          mat1(io1,io2,ispin,ik) = cbox * system%hvol
        end do
        end do
        end do
        end do
      end if

      call comm_summation(mat1,mat,no0*no*nspin*nk,info%icomm_rko)

    end subroutine inner_product
    
    subroutine write_intraband_current
      use density_matrix, only:calc_current_decomposed
      use inputoutput, only: t_unit_current
      implicit none
      real(8) :: curr_intra(3,nspin)
      real(8) :: curr_decomp(3,nspin,no0,nk)
      
      call calc_current_decomposed(rt%system_proj,mg,stencil,rt%info_proj,rt%srg_proj,rt%tpsi0,ppg,curr_decomp)
      
      curr_intra = 0d0
      do ispin=1,nspin_tmp
        do ik=1,nk
          do io1=1,no0
            do io2=1,no
              curr_intra(1:3,ispin) = curr_intra(1:3,ispin) + system%rocc(io2,ik,ispin)*system%wtk(ik) &
              & * abs(mat(io1,io2,ispin,ik))**2 * curr_decomp(1:3,ispin,io1,ik)
            end do
          end do
        end do
      end do
      
      if(comm_is_root(nproc_id_global))then
        write(ofl%fh_intra_current, "(F16.8)",advance='no') itt * dt * t_unit_time%conv
        if(nspin_tmp==1) then
          write(ofl%fh_intra_current, "(99(1X,E23.15E3))") &
          & curr_intra(1:3,1) * t_unit_current%conv
        else if(nspin_tmp==2) then
          write(ofl%fh_intra_current, "(99(1X,E23.15E3))") &
          & curr_intra(1:3,1) * t_unit_current%conv, &
          & curr_intra(1:3,2) * t_unit_current%conv
        end if
      end if
    end subroutine write_intraband_current
    
  end subroutine projection

!===================================================================================================================================

  subroutine write_current_decomposed(itt,ofl,mg,system,info,stencil,srg,psi,ppg)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: yn_spinorbit,base_directory,SYSname
    use inputoutput, only: t_unit_current
    use filesystem, only: open_filehandle
    use density_matrix, only:calc_current_decomposed
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)                       :: ofl
    type(s_rgrid)           ,intent(in) :: mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_orbital)                     :: psi
    type(s_sendrecv_grid)               :: srg
    integer :: ik,io,ispin,nspin_tmp
    character(256) :: file_current_decomposed
    real(8) :: curr_decomp(3,system%nspin,system%no,system%nk)
    
    !(header in SYSname_current_decomposed.data)
    if(itt < 0) then
      if(comm_is_root(nproc_id_global))then
        write(file_current_decomposed,"(2A,'_current_decomposed.data')") trim(base_directory),trim(SYSname)
        ofl%fh_current_decomposed = open_filehandle(file_current_decomposed)
        open(ofl%fh_current_decomposed,file=file_current_decomposed)
        write(ofl%fh_current_decomposed, '("#",1X,A)') "decomposition of the current density"
        write(ofl%fh_current_decomposed, '("#",1X,A,":",1X,A)') "it", "time step index"
        write(ofl%fh_current_decomposed, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(ofl%fh_current_decomposed, '("#",1X,A,":",1X,A)') "ib", "band index"
        write(ofl%fh_current_decomposed, '("#",1X,A,":",1X,A)') "is", "spin index"
        write(ofl%fh_current_decomposed, '("#",1X,A,":",1X,A)') "curr_decomp_[xyz]", "decomposed elements the current density"
        write(ofl%fh_current_decomposed, '("#",1X,A)') "------------------------------------------------"
        write(ofl%fh_current_decomposed, '("#",1X,A)') "it"
        write(ofl%fh_current_decomposed, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "ib", "none", &
          & 3, "is", "none", &
          & 4, "curr_decomp_x", trim(t_unit_current%name), &
          & 5, "curr_decomp_y", trim(t_unit_current%name), &
          & 6, "curr_decomp_z", trim(t_unit_current%name)
      end if
      return
    end if
    
    nspin_tmp = system%nspin
    if(yn_spinorbit=='y') then
      nspin_tmp = 1
    end if
    
    call calc_current_decomposed(system,mg,stencil,info,srg,psi,ppg,curr_decomp)
    
    if(comm_is_root(nproc_id_global))then
      write(ofl%fh_current_decomposed,'(i11)') itt
      do ik=1,system%nk
      do io=1,system%no
      do ispin=1,nspin_tmp
        write(ofl%fh_current_decomposed,'(3(1X,i6),3(1X,E23.15E3))') &
        & ik,io,ispin, curr_decomp(1:3,ispin,io,ik)* (system%rocc(io,ik,ispin)*system%wtk(ik)) * t_unit_current%conv
      end do
      end do
      end do
    end if
    
  end subroutine write_current_decomposed

!===================================================================================================================================

  subroutine write_rt_spin(itt,ofl,system,lg,mg,info,stencil,ppg,psi)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: dt,yn_spinorbit,base_directory,SYSname, &
    & yn_out_mag_decomposed_rt, yn_out_mag_micro_rt, &
    & yn_out_spin_current_decomposed,yn_out_spin_current_micro
    use filesystem, only: open_filehandle
    use inputoutput, only: t_unit_current,t_unit_time
    use noncollinear_module, only: calc_magnetization,calc_magnetization_decomposed,calc_spin_current,calc_magnetization_micro
    use writefield, only: write_spin_current_micro,write_magnetization_micro
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)                       :: ofl
    type(s_dft_system)      ,intent(in) :: system
    type(s_rgrid)           ,intent(in) :: lg,mg
    type(s_parallel_info)   ,intent(in) :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_orbital)         ,intent(in) :: psi
    !
    integer ik,io,i
    real(8) :: m(3),mag_orb(3,system%no,system%nk),spin_curr(3,0:3)
    real(8) :: m_micro(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:3)
    real(8) :: spin_curr_micro(3,0:3, &
                & mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    real(8) :: spin_curr_band(3,0:3,system%no,system%nk)
    
    if(yn_spinorbit=='n') stop "error: write_rt_spin with yn_spinorbit=n"
    
    if(itt < 0) then
      call file_header
      return
    end if
    
    call calc_magnetization(system,mg,info,m)
    if(yn_out_mag_decomposed_rt=='y') call calc_magnetization_decomposed(system,mg,info,psi,mag_orb)
    if(yn_out_mag_micro_rt=='y') call calc_magnetization_micro(mg,m_micro)

    call calc_spin_current(system,mg,stencil,info,psi,ppg,spin_curr_micro,spin_curr_band)
    
    spin_curr = 0d0
    do ik=1,system%nk
    do io=1,system%no
      spin_curr = spin_curr + spin_curr_band(:,:,io,ik)* system%rocc(io,ik,1)*system%wtk(ik)
    end do
    end do
  
    if(comm_is_root(nproc_id_global))then
    
      write(ofl%fh_rt_spin,'(1000(1X,E23.15E3))') dble(itt)*dt*t_unit_time%conv,m(1),m(2),m(3), &
      & spin_curr(1:3,0) * t_unit_current%conv, &
      & spin_curr(1:3,1) * t_unit_current%conv, &
      & spin_curr(1:3,2) * t_unit_current%conv, &
      & spin_curr(1:3,3) * t_unit_current%conv
      
      if(yn_out_mag_decomposed_rt=='y') then
        write(ofl%fh_mag_decomposed_rt,'(i11)') itt
        do ik=1,system%nk
        do io=1,system%no
          write(ofl%fh_mag_decomposed_rt,'(i6,1X,i6,1000(1X,E23.15E3))') ik,io,mag_orb(1,io,ik),mag_orb(2,io,ik),mag_orb(3,io,ik)
        end do
        end do
      end if

      if(yn_out_spin_current_decomposed=='y') then
        write(ofl%fh_spin_current_decomposed,'(i11)') itt
        do ik=1,system%nk
        do io=1,system%no
        do i=0,3
          write(ofl%fh_spin_current_decomposed,'(3(1X,i6),3(1X,E23.15E3))') &
          & ik,io,i, spin_curr_band(1:3,i,io,ik)* (system%rocc(io,ik,1)*system%wtk(ik)) * t_unit_current%conv
        end do
        end do
        end do
      end if
      
    end if
  
    if(yn_out_mag_micro_rt=='y') then
      call write_magnetization_micro(lg,mg,system,info,itt,m_micro)
    end if
    
    if(yn_out_spin_current_micro=='y') then
      call write_spin_current_micro(lg,mg,system,info,itt,spin_curr_micro)
    end if
    
  contains
  
    subroutine file_header
      implicit none
      integer :: fh
      character(256) :: filename
      !
      if(comm_is_root(nproc_id_global))then
        !(header in SYSname_rt_spin.data)
        write(filename,"(2A,'_rt_spin.data')") trim(base_directory),trim(SYSname)
        ofl%fh_rt_spin = open_filehandle(filename)
        fh = ofl%fh_rt_spin
        open(fh,file=filename)
        write(fh, '("#",1X,A)') "Spin magnetization & Spin current density"
        write(fh, '("#",1X,A,":",1X,A)') "mag_[xyz]", "Total magnetization"
        write(fh, '("#",1X,A,":",1X,A)') "spin_curr_i_[xyz]", "Spin current density with the Pauli matrix sigma_i (i=0,x,y,z)"
        write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "Time", trim(t_unit_time%name), &
        & 2,  "mag_x", "none", &
        & 3,  "mag_y", "none", &
        & 4,  "mag_z", "none", &
        & 5,  "spin_curr_0_x", trim(t_unit_current%name), &
        & 6,  "spin_curr_0_y", trim(t_unit_current%name), &
        & 7,  "spin_curr_0_z", trim(t_unit_current%name), &
        & 8,  "spin_curr_x_x", trim(t_unit_current%name), &
        & 9,  "spin_curr_x_y", trim(t_unit_current%name), &
        & 10, "spin_curr_x_z", trim(t_unit_current%name), &
        & 11, "spin_curr_y_x", trim(t_unit_current%name), &
        & 12, "spin_curr_y_y", trim(t_unit_current%name), &
        & 13, "spin_curr_y_z", trim(t_unit_current%name), &
        & 14, "spin_curr_z_x", trim(t_unit_current%name), &
        & 15, "spin_curr_z_y", trim(t_unit_current%name), &
        & 16, "spin_curr_z_z", trim(t_unit_current%name)
          
        !(header in SYSname_mag_decomposed_rt.data)
        if(yn_out_mag_decomposed_rt=='y') then
          write(filename,"(2A,'_mag_decomposed_rt.data')") trim(base_directory),trim(SYSname)
          ofl%fh_mag_decomposed_rt = open_filehandle(filename)
          fh = ofl%fh_mag_decomposed_rt
          open(fh,file=filename)
          write(fh, '("#",1X,A)') "band decomposition of the time-dependent spin magnetization"
          write(fh, '("#",1X,A,":",1X,A)') "it", "time step index"
          write(fh, '("#",1X,A,":",1X,A)') "ik", "k-point index"
          write(fh, '("#",1X,A,":",1X,A)') "ib", "band index"
          write(fh, '("#",1X,A,":",1X,A)') "mag_orb_[xyz]", "decomposed elements of the magnetization"
          write(fh, '("#",1X,A)') "------------------------------------------------"
          write(fh, '("#",1X,A)') "it"
          write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "ib", "none", &
          & 3, "mag_orb_x", "none", &
          & 4, "mag_orb_y", "none", &
          & 5, "mag_orb_z", "none"
        end if
        
        !(header in SYSname_spin_current_decomposed.data)
        if(yn_out_spin_current_decomposed=='y') then
          write(filename,"(2A,'_spin_current_decomposed.data')") trim(base_directory),trim(SYSname)
          ofl%fh_spin_current_decomposed = open_filehandle(filename)
          fh = ofl%fh_spin_current_decomposed
          open(fh,file=filename)
          write(fh, '("#",1X,A)') "Band decomposition of the spin current density"
          write(fh, '("#",1X,A,":",1X,A)') "it", "time step index"
          write(fh, '("#",1X,A,":",1X,A)') "ik", "k-point index"
          write(fh, '("#",1X,A,":",1X,A)') "ib", "band index"
          write(fh, '("#",1X,A,":",1X,A)') "i", "index of the Pauli matrices (0,x,y,z)"
          write(fh, '("#",1X,A,":",1X,A)') "spin_curr_[xyz]", "decomposed elements of the spin current density"
          write(fh, '("#",1X,A)') "------------------------------------------------"
          write(fh, '("#",1X,A)') "it"
          write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "ib", "none", &
          & 3, "i", "none", &
          & 4, "spin_curr_x", trim(t_unit_current%name), &
          & 5, "spin_curr_y", trim(t_unit_current%name), &
          & 6, "spin_curr_z", trim(t_unit_current%name)
        end if
        
      end if
    
    end subroutine file_header
    
  end subroutine write_rt_spin
  
  subroutine write_mag_decomposed_gs(system,mg,info,psi)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: yn_spinorbit, base_directory, sysname
    use noncollinear_module, only: calc_magnetization,calc_magnetization_decomposed
    use filesystem, only: open_filehandle
    implicit none
    type(s_dft_system)      ,intent(in) :: system
    type(s_rgrid)           ,intent(in) :: mg
    type(s_parallel_info)   ,intent(in) :: info
    type(s_orbital)         ,intent(in) :: psi
    !
    integer ik,io,fh_gs_mag
    character(256) :: file_gs_mag
    real(8) :: m(3),mag_orb(3,system%no,system%nk)
    
    if(yn_spinorbit=='n') stop "error: write_mag_decomposed_gs with yn_spinorbit=n"
    
    call calc_magnetization(system,mg,info,m)
    call calc_magnetization_decomposed(system,mg,info,psi,mag_orb)
    
    if(comm_is_root(nproc_id_global))then
      !(header in gs_mag.data)
      write(file_gs_mag,"(2A,'_mag_decomposed_gs.data')") trim(base_directory),trim(SYSname)
      fh_gs_mag = open_filehandle(file_gs_mag)
      open(fh_gs_mag,file=file_gs_mag)
      write(fh_gs_mag, '("#",1X,A)') "Band decomposition of the spin magnetization @ the ground state"
      write(fh_gs_mag, '("#",1X,A,":",1X,A)') "ik", "k-point index"
      write(fh_gs_mag, '("#",1X,A,":",1X,A)') "ib", "band index"
      write(fh_gs_mag, '("#",1X,A,":",1X,A)') "mag", "Total magnetization"
      write(fh_gs_mag, '("#",1X,A,":",1X,A)') "mag_orb", "decomposed elements of the magnetization"
      write(fh_gs_mag, '("#",99(1X,I0,":",A,"[",A,"]"))') &
      & 1, "mag(1)", "none", &
      & 2, "mag(2)", "none", &
      & 3, "mag(3)", "none"
      write(fh_gs_mag, '("#",99(1X,I0,":",A,"[",A,"]"))') &
      & 1, "ik", "none", &
      & 2, "io", "none", &
      & 3, "mag_orb(1)", "none", &
      & 4, "mag_orb(2)", "none", &
      & 5, "mag_orb(3)", "none"
      write(fh_gs_mag,'(1000(1X,E23.15E3))') m(1),m(2),m(3)
      do ik=1,system%nk
      do io=1,system%no
        write(fh_gs_mag,'(i6,1X,i6,1000(1X,E23.15E3))') ik,io,mag_orb(1,io,ik),mag_orb(2,io,ik),mag_orb(3,io,ik)
      end do
      end do
      close(fh_gs_mag)
    end if
    
  end subroutine write_mag_decomposed_gs

!===================================================================================================================================

end module write_sub
