!
!  Copyright 2019 SALMON developers
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
!-----------------------------------------------------------------------------------------
subroutine eh_init(fs,fw)
  use inputoutput,          only: nt_em,al_em,dl_em,dt_em,iboundary,&
                                  utime_from_au,ulength_from_au,uenergy_from_au,unit_system,&
                                  uenergy_to_au,ulength_to_au,ucharge_to_au,iperiodic,directory,&
                                  imedia_num,shape_file,epsilon,rmu,sigma,type_media,&
                                  omega_p_d,gamma_d,smooth_d,weight_d,&
                                  iobs_num_em,obs_loc_em,wave_input,trans_longi,e_impulse,nenergy,&
                                  source_loc1,ek_dir1,epdir_re1,epdir_im1,ae_shape1,&
                                  phi_cep1,rlaser_int_wcm2_1,amplitude1,&
                                  source_loc2,ek_dir2,epdir_re2,epdir_im2,ae_shape2,&
                                  phi_cep2,rlaser_int_wcm2_2,amplitude2
  use salmon_parallel,      only: nproc_id_global, nproc_group_global
  use salmon_communication, only: comm_is_root, comm_bcast
  use structures,           only: s_fdtd_system
  use salmon_maxwell,       only: s_fdtd_work
  implicit none
  type(s_fdtd_system) :: fs
  type(s_fdtd_work)   :: fw
  integer             :: ii,ij,ix,iy,iz,icount,icount_d,iflag
  real(8),parameter   :: pi=3.141592653589793d0
  real(8)             :: dt_cfl,diff_cep
  character(1)        :: dir
  character(128)      :: save_name
  
  !set initial parameter and value
  fw%c_0        = 1.370359991378353d2
  fw%Nd         = 1
  fw%iter_sta   = 1
  fw%iter_end   = nt_em
  fs%rlsize(:) = al_em(:)
  fs%hgs(:)    = dl_em(:)
  fw%ifn        = 600
  fw%ipml_l     = 8
  fw%pml_m      = 4.0d0
  fw%pml_r      = 1.0d-7
  if(iperiodic==0) then
    fs%i_bc(:,:)=1 !PML
    do ix=1,3
    do iy=1,2
      if(iboundary(ix,iy)==1) then
        fs%i_bc(ix,iy)=0 !PEC
      end if
    end do
    end do
  elseif(iperiodic==3) then
    fs%i_bc(:,:)  = iboundary(:,:) !Periodic or PML
  end if
  select case(unit_system)
  case('au','a.u.')
    fw%uVperm_from_au=1.0d0
    fw%uAperm_from_au=1.0d0
  case('A_eV_fs')
    !see amplitude1 or amplitude2 in src/io/iunputoutput.f90
    fw%uVperm_from_au=1/(uenergy_to_au/ulength_to_au/ucharge_to_au)
    fw%uAperm_from_au=fw%uVperm_from_au
  end select
  
  !prepare GCEED(set mpi condition, gird, and sendrecv environment)
  call eh_prep_GCEED(fs,fw)
  
  !set coordinate
  allocate(fw%coo(minval(fs%lg_sta(:))-fw%Nd:maxval(fs%lg_end(:))+fw%Nd,3))
  call eh_set_coo(iperiodic,fw%Nd,fw%ioddeven(:),fs%lg_sta(:),fs%lg_end(:),fs%hgs(:),fw%coo)
  
  !set and check dt
  dt_cfl=1.0d0/( &
         fw%c_0*sqrt( (1.0d0/fs%hgs(1))**2.0d0+(1.0d0/fs%hgs(2))**2.0d0+(1.0d0/fs%hgs(3))**2.0d0 ) &
         )
  if(dt_em==0.0d0) then
    fs%dt=dt_cfl*0.99d0
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "**************************"
      write(*,*) "From CFL condition, dt_em is determined by", fs%dt*utime_from_au
      write(*,*) "in the unit system, ",trim(unit_system),"."
      write(*,*) "**************************"
    end if
  elseif(dt_em>=dt_cfl) then
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "**************************"
      write(*,*) "To sufficient CFL condition, dt_em must be set smaller than", dt_cfl*utime_from_au
      write(*,*) "in the unit system, ",trim(unit_system),"."
      write(*,*) "**************************"
    end if
    stop
  else
    fs%dt=dt_em
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "**************************"
      write(*,*) "dt_em =", fs%dt*utime_from_au
      write(*,*) "in the unit system, ",trim(unit_system),"."
      write(*,*) "**************************"
    end if
  end if
  call comm_bcast(fs%dt,nproc_group_global)
  
  !allocate
  call eh_allocate
  allocate(fw%c2_jx(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                    fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                    fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
           fw%c2_jy(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                    fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                    fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
           fw%c2_jz(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                    fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                    fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd) )
  fw%c2_jx(:,:,:)=0.0d0; fw%c2_jy(:,:,:)=0.0d0; fw%c2_jz(:,:,:)=0.0d0;
  
  !input fdtd shape
  allocate(fs%imedia(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
  fs%imedia(:,:,:)=0
  allocate(fw%rmedia(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
  fw%rmedia(:,:,:)=0.0d0
  if(imedia_num>0) then
    !check file format and input shape file
    if(comm_is_root(nproc_id_global)) write(*,*)
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
    if(index(shape_file,".cube", back=.true.)/=0) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file is inputed by .cube format."
      end if
      call eh_input_shape(fw%ifn,fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fs%imedia,'cu')
      fw%rmedia(:,:,:)=dble(fs%imedia(:,:,:))
      call eh_sendrecv(fs,fw,'r')
      fs%imedia(:,:,:)=int(fw%rmedia(:,:,:)+1d-3)
    elseif(index(shape_file,".mp", back=.true.)/=0) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file is inputed by .mp format."
        write(*,*) "This version works for only .cube format.."
      end if
      stop
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file must be .cube or .mp formats."
      end if
      stop
    end if
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
  end if
  
  !prepare Drude
  fw%inum_d=0
  do ii=0,imedia_num
    select case(type_media(ii))
    case('DRUDE','Drude','drude','D','d')
      fw%inum_d=fw%inum_d+1
    end select
  end do
  if(fw%inum_d>0) then
    !set counter and make imedia_d
    icount_d=1
    allocate(fw%imedia_d(fw%inum_d))
    fw%imedia_d(:)=0;
    do ii=0,imedia_num
      select case(type_media(ii))
      case('DRUDE','Drude','drude','D','d')
        fw%imedia_d(icount_d)=ii
        icount_d=icount_d+1;
      end select
    end do
    
    !reset counter
    icount_d=1
    
    !allocate drude variable
    allocate(fw%idx_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                      fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                      fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d),&
             fw%idy_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                      fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                      fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d),&
             fw%idz_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                      fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                      fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d) )
    fw%idx_d(:,:,:,:)=0; fw%idy_d(:,:,:,:)=0; fw%idz_d(:,:,:,:)=0;
    allocate( fw%rjx_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d),&
              fw%rjy_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d),&
              fw%rjz_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d) )
    fw%rjx_d(:,:,:,:)=0.0d0; fw%rjy_d(:,:,:,:)=0.0d0; fw%rjz_d(:,:,:,:)=0.0d0;
    allocate( fw%rjx_sum_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                           fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                           fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
              fw%rjy_sum_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                           fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                           fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
              fw%rjz_sum_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                           fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                           fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd) )
    fw%rjx_sum_d(:,:,:)=0.0d0; fw%rjy_sum_d(:,:,:)=0.0d0; fw%rjz_sum_d(:,:,:)=0.0d0;
    allocate( fw%wex_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d),&
              fw%wey_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d),&
              fw%wez_d(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd,fw%inum_d) )
    fw%wex_d(:,:,:,:)=0.0d0; fw%wey_d(:,:,:,:)=0.0d0; fw%wez_d(:,:,:,:)=0.0d0;
    allocate(fw%c1_j_d(fw%inum_d),fw%c2_j_d(fw%inum_d))
    fw%c1_j_d(:)=0.0d0; fw%c2_j_d(:)=0.0d0;
  end if
  
  !set fdtd coeffient and write media information
  allocate(fw%rep(0:imedia_num),fw%rmu(0:imedia_num),fw%sig(0:imedia_num))
  fw%rep(:)=1.0d0; fw%rmu(:)=1.0d0; fw%sig(:)=0.0d0;
  do ii=0,imedia_num
    call eh_coeff
  end do
  deallocate(fs%imedia); deallocate(fw%rmedia);
  if(comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "**************************"
    write(*,'(A,I3)')          ' imedia_num = ',imedia_num
    do ii=0,imedia_num
      write(*,*) "=========================="
      write(*,'(A,I3,A)')      ' id ',ii, ':'
      select case(type_media(ii))
      case ('VACUUM','Vacuum','vacuum')
        if(epsilon(ii)/=1d0 .or. rmu(ii)/=1d0 .or. sigma(ii)/=0d0) then
          write(*,'(A)'  )     ' type_media =  constant media'
        else
          write(*,'(A,A)')       ' type_media =  ', trim(type_media(ii))
        end if
      case('DRUDE','Drude','drude','D','d')
        write(*,'(A,A)')       ' type_media =  ', trim(type_media(ii))
        write(*,'(A,ES12.5)')  ' omega_p_d  = ', omega_p_d(ii)*uenergy_from_au
        write(*,'(A,ES12.5)')  ' gamma_d    = ', gamma_d(ii)*uenergy_from_au
      case default
        write(*,'(A,A)')       ' type_media =  ', trim(type_media(ii))
      end select
      write(*,'(A,ES12.5)')    ' epsilon    = ', epsilon(ii)
      write(*,'(A,ES12.5)')    ' rmu        = ', rmu(ii)
      write(*,'(A,ES12.5)'   ) ' sigma      = ', sigma(ii)
    end do
    write(*,*) "**************************"
  end if
  
  !apply smoothing
  call eh_smoothing
  
  !set calculation area
  fw%iex_y_sta(:)=fs%ng_sta(:); fw%iex_y_end(:)=fs%ng_end(:);
  fw%iex_z_sta(:)=fs%ng_sta(:); fw%iex_z_end(:)=fs%ng_end(:);
  fw%iey_z_sta(:)=fs%ng_sta(:); fw%iey_z_end(:)=fs%ng_end(:);
  fw%iey_x_sta(:)=fs%ng_sta(:); fw%iey_x_end(:)=fs%ng_end(:);
  fw%iez_x_sta(:)=fs%ng_sta(:); fw%iez_x_end(:)=fs%ng_end(:);
  fw%iez_y_sta(:)=fs%ng_sta(:); fw%iez_y_end(:)=fs%ng_end(:);
  fw%ihx_y_sta(:)=fs%ng_sta(:); fw%ihx_y_end(:)=fs%ng_end(:);
  fw%ihx_z_sta(:)=fs%ng_sta(:); fw%ihx_z_end(:)=fs%ng_end(:);
  fw%ihy_z_sta(:)=fs%ng_sta(:); fw%ihy_z_end(:)=fs%ng_end(:);
  fw%ihy_x_sta(:)=fs%ng_sta(:); fw%ihy_x_end(:)=fs%ng_end(:);
  fw%ihz_x_sta(:)=fs%ng_sta(:); fw%ihz_x_end(:)=fs%ng_end(:);
  fw%ihz_y_sta(:)=fs%ng_sta(:); fw%ihz_y_end(:)=fs%ng_end(:);
  if((fs%i_bc(1,1)==1).and.(fs%ng_sta(1)==fs%lg_sta(1))) then !x, bottom
    fw%iey_x_sta(1)=fs%ng_sta(1)+1; fw%iez_x_sta(1)=fs%ng_sta(1)+1;
  end if
  if((fs%i_bc(1,2)==1).and.(fs%ng_end(1)==fs%lg_end(1))) then !x, top
    fw%iex_y_end(1)=fs%ng_end(1)-1; fw%iex_z_end(1)=fs%ng_end(1)-1;
    fw%iey_x_end(1)=fs%ng_end(1)-1; fw%iez_x_end(1)=fs%ng_end(1)-1;
    fw%ihy_z_end(1)=fs%ng_end(1)-1; fw%ihy_x_end(1)=fs%ng_end(1)-1;
    fw%ihz_x_end(1)=fs%ng_end(1)-1; fw%ihz_y_end(1)=fs%ng_end(1)-1;
  end if
  if((fs%i_bc(2,1)==1).and.(fs%ng_sta(2)==fs%lg_sta(2))) then !y, bottom
    fw%iex_y_sta(2)=fs%ng_sta(2)+1; fw%iez_y_sta(2)=fs%ng_sta(2)+1;
  end if
  if((fs%i_bc(2,2)==1).and.(fs%ng_end(2)==fs%lg_end(2))) then !y, top
    fw%iex_y_end(2)=fs%ng_end(2)-1; fw%iey_z_end(2)=fs%ng_end(2)-1;
    fw%iey_x_end(2)=fs%ng_end(2)-1; fw%iez_y_end(2)=fs%ng_end(2)-1;
    fw%ihx_y_end(2)=fs%ng_end(2)-1; fw%ihx_z_end(2)=fs%ng_end(2)-1;
    fw%ihz_x_end(2)=fs%ng_end(2)-1; fw%ihz_y_end(2)=fs%ng_end(2)-1;
  end if
  if((fs%i_bc(3,1)==1).and.(fs%ng_sta(3)==fs%lg_sta(3))) then !z, bottom
    fw%iex_z_sta(3)=fs%ng_sta(3)+1; fw%iey_z_sta(3)=fs%ng_sta(3)+1;
  end if
  if((fs%i_bc(3,2)==1).and.(fs%ng_end(3)==fs%lg_end(3))) then !z, top
    fw%iex_z_end(3)=fs%ng_end(3)-1; fw%iey_z_end(3)=fs%ng_end(3)-1;
    fw%iez_x_end(3)=fs%ng_end(3)-1; fw%iez_y_end(3)=fs%ng_end(3)-1;
    fw%ihx_y_end(3)=fs%ng_end(3)-1; fw%ihx_z_end(3)=fs%ng_end(3)-1;
    fw%ihy_z_end(3)=fs%ng_end(3)-1; fw%ihy_x_end(3)=fs%ng_end(3)-1;
  end if
  
  !set pml
  call eh_set_pml(1,fw%c1_ey_x,fw%c2_ey_x,fw%c1_ez_x,fw%c2_ez_x,&
                    fw%c1_hy_x,fw%c2_hy_x,fw%c1_hz_x,fw%c2_hz_x) !x direction
  call eh_set_pml(2,fw%c1_ez_y,fw%c2_ez_y,fw%c1_ex_y,fw%c2_ex_y,&
                    fw%c1_hz_y,fw%c2_hz_y,fw%c1_hx_y,fw%c2_hx_y) !y direction
  call eh_set_pml(3,fw%c1_ex_z,fw%c2_ex_z,fw%c1_ey_z,fw%c2_ey_z,&
                    fw%c1_hx_z,fw%c2_hx_z,fw%c1_hy_z,fw%c2_hy_z) !z direction
  if(maxval(fs%i_bc(:,:))>0) then
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      do ii=1,3
        if(ii==1) then
          dir='x'
        elseif(ii==2) then
          dir='y'
        elseif(ii==3) then
          dir='z'
        end if
        if(fs%i_bc(ii,1)==1) write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                               ' PML has been set for ',dir,'-direction: ',&
                               fw%coo(fs%lg_sta(ii),ii)*ulength_from_au,' to ',&
                               fw%coo(fs%lg_sta(ii)+fw%ipml_l,ii)*ulength_from_au,'.'
        if(fs%i_bc(ii,2)==1) write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                               ' PML has been set for ',dir,'-direction: ',&
                               fw%coo(fs%lg_end(ii)-fw%ipml_l,ii)*ulength_from_au,' to ',&
                               fw%coo(fs%lg_end(ii),ii)*ulength_from_au,'.'
      end do
      write(*,*) "**************************"
    end if
  end if
  
  !prepare observation
  if(iobs_num_em>0) then
    !set initial
    allocate(fw%ex_s(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%ey_s(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%ez_s(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%hx_s(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%hy_s(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%hz_s(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ex_s(:,:,:)=0; fw%ey_s(:,:,:)=0; fw%ez_s(:,:,:)=0; 
    fw%hx_s(:,:,:)=0; fw%hy_s(:,:,:)=0; fw%hz_s(:,:,:)=0; 
    allocate(fw%iobs_po_id(iobs_num_em,3)) !1:x,        2:y,        3:z
    allocate(fw%iobs_po_pe(iobs_num_em))
    allocate(fw%iobs_li_pe(iobs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
    allocate(fw%iobs_pl_pe(iobs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
    fw%iobs_po_id(:,:)=0; fw%iobs_po_pe(:)=0; fw%iobs_li_pe(:,:)=0; fw%iobs_pl_pe(:,:)=0; 
    fw%e_max=0.0d0; fw%h_max=0.0d0;
    
    !search observation point
    do ii=1,iobs_num_em
      call eh_find_point(obs_loc_em(ii,:),fw%iobs_po_id(ii,:),&
                         fw%iobs_po_pe(ii),fw%iobs_li_pe(ii,:),fw%iobs_pl_pe(ii,:),fs%ng_sta,fs%ng_end,&
                         minval(fs%lg_sta)-fw%Nd,maxval(fs%lg_end)+fw%Nd,fw%coo)
    end do
    
    !write information
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      if(iobs_num_em==1) then
        write(*,*) "Observation point is placed at"
      else
        write(*,*) "Observation points are placed at"
      end if
      do ii=1,iobs_num_em
        write(*,'(I3,A,3ES14.5)') ii,":",(fw%coo(fw%iobs_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
      end do
      write(*,*) "**************************"
      do ii=1,iobs_num_em
        write(save_name,*) ii
        save_name=trim(adjustl(directory))//'/obs'//trim(adjustl(save_name))//'_at_point.data'
        open(fw%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'(A)') "# time[a.u.], Ex[a.u.], Ey[a.u.], Ez[a.u.], Hx[a.u.], Hy[a.u.], Hz[a.u.]" 
        case('A_eV_fs')
          write(fw%ifn,'(A)') "# time[fs], Ex[V/Ang.], Ey[V/Ang.], Ez[V/Ang.], Hx[A/Ang.], Hy[A/Ang.], Hz[A/Ang.]" 
        end select
        close(fw%ifn)
      end do
    end if
  end if
  
  !check incident current source condition
  select case(wave_input)
  case('source')
    !linear response
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ae_shape1/2:"
        write(*,*) "For ae_shape1/2 = impulse, wave_input must be default(do not set source)."
      end if
      stop
    end if
    
    !source1
    if    (ek_dir1(1)==0.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==0.0d0) then 
      fw%inc_dist1='none'
    elseif(ek_dir1(1)==1.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==0.0d0) then
      if(epdir_re1(1)/=0.0d0.or.epdir_im1(1)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(1) and epdir_im1(1):"
          write(*,*) "For theory = Maxwell and ek_dir1(1) = 1.0d0, epdir_re1(1) and epdir_im1(1) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist1='yz-plane'
      end if
    elseif(ek_dir1(1)==0.0d0.and.ek_dir1(2)==1.0d0.and.ek_dir1(3)==0.0d0) then
      if(epdir_re1(2)/=0.0d0.or.epdir_im1(2)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(2) and epdir_im1(2):"
          write(*,*) "For theory = Maxwell and ek_dir1(2) = 1.0d0, epdir_re1(2) and epdir_im1(2) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist1='xz-plane'
      end if
    elseif(ek_dir1(1)==0.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==1.0d0) then
      if(epdir_re1(3)/=0.0d0.or.epdir_im1(3)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(3) and epdir_im1(3):"
          write(*,*) "For theory = Maxwell and ek_dir1(3) = 1.0d0, epdir_re1(3) and epdir_im1(3) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist1='xy-plane'
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ek_dir1:"
        write(*,*) "For theory = Maxwell, ek_dir1 is only allowed by"
        write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
      end if
      stop
    end if
    
    !source2
    if    (ek_dir2(1)==0.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==0.0d0) then 
      fw%inc_dist2='none'
    elseif(ek_dir2(1)==1.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==0.0d0) then
      if(epdir_re2(1)/=0.0d0.or.epdir_im2(1)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(1) and epdir_im2(1):"
          write(*,*) "For theory = Maxwell and ek_dir2(1) = 1.0d0, epdir_re2(1) and epdir_im2(1) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist2='yz-plane'
      end if
    elseif(ek_dir2(1)==0.0d0.and.ek_dir2(2)==1.0d0.and.ek_dir2(3)==0.0d0) then
      if(epdir_re2(2)/=0.0d0.or.epdir_im2(2)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(2) and epdir_im2(2):"
          write(*,*) "For theory = Maxwell and ek_dir2(2) = 1.0d0, epdir_re2(2) and epdir_im2(2) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist2='xz-plane'
      end if
    elseif(ek_dir2(1)==0.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==1.0d0) then
      if(epdir_re2(3)/=0.0d0.or.epdir_im2(3)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(3) and epdir_im2(3):"
          write(*,*) "For theory = Maxwell and ek_dir2(3) = 1.0d0, epdir_re2(3) and epdir_im2(3) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist2='xy-plane'
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ek_dir2:"
        write(*,*) "For theory = Maxwell, ek_dir2 is only allowed by"
        write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
      end if
      stop
    end if
  case('point','x-line','y-line','z-line')
    !these selection are for debug
    fw%inc_dist1=wave_input; fw%inc_dist2='none';
    if(comm_is_root(nproc_id_global)) write(*,*) trim(wave_input), " source is used."
  case default
    fw%inc_dist1='none'; fw%inc_dist2='none';
    if(ae_shape1/='impulse'.and.ae_shape2/='impulse') then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid wave_input:"
        write(*,*) "For theory = Maxwell, wave_input must be source"
        write(*,*) "or ae_shape1 and/or ae_shape2 must be impulse."
      end if
      stop
    end if
  end select
  
  !prepare incident current source
  if((fw%inc_dist1=='none').and.(fw%inc_dist2=='none')) then
    fw%inc_num=0
  else
    fw%inc_num=2
  end if
  if(fw%inc_num>0) then
    !set initial
    allocate(fw%inc_po_id(iobs_num_em,3)) !1:x,        2:y,        3:z
    allocate(fw%inc_po_pe(iobs_num_em))
    allocate(fw%inc_li_pe(iobs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
    allocate(fw%inc_pl_pe(iobs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
    fw%inc_po_id(:,:)=0; fw%inc_po_pe(:)=0; fw%inc_li_pe(:,:)=0; fw%inc_pl_pe(:,:)=0; 
    do ii=1,3
      fw%c2_inc_xyz(ii)=(fw%c_0/fw%rep(0)*fs%dt) &
                         /(1.0d0+2.0d0*pi*fw%sig(0)/fw%rep(0)*fs%dt) &
                         *2.0d0/( fs%hgs(ii)*sqrt(fw%rmu(0)/fw%rep(0)) )
    end do
    
    !search incident current source point and check others
    if(fw%inc_dist1/='none') then
      ii=1
      call eh_find_point(source_loc1(:),fw%inc_po_id(ii,:),&
                         fw%inc_po_pe(ii),fw%inc_li_pe(ii,:),fw%inc_pl_pe(ii,:),fs%ng_sta,fs%ng_end,&
                         minval(fs%lg_sta(:))-fw%Nd,maxval(fs%lg_end(:))+fw%Nd,fw%coo(:,:))
      select case(ae_shape1)
      case("Ecos2","Acos2")
        continue
      case default
        stop 'set ae_shape1 to "Ecos2" or "Acos2"'
      end select
      diff_cep=(phi_cep1-0.25d0)*2.d0-int((phi_cep1-0.25d0)*2.d0)
      if(ae_shape1=="Ecos2".and.abs(diff_cep)>=1.d-12)then
        stop "phi_cep1 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape1."
      end if
      if(rlaser_int_wcm2_1/=-1d0) &
        amplitude1=sqrt(rlaser_int_wcm2_1)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    if(fw%inc_dist2/='none') then
      ii=2
      call eh_find_point(source_loc2(:),fw%inc_po_id(ii,:),&
                         fw%inc_po_pe(ii),fw%inc_li_pe(ii,:),fw%inc_pl_pe(ii,:),fs%ng_sta,fs%ng_end,&
                         minval(fs%lg_sta(:))-fw%Nd,maxval(fs%lg_end(:))+fw%Nd,fw%coo(:,:))
      select case(ae_shape2)
      case("Ecos2","Acos2")
        continue
      case default
        stop 'set ae_shape2 to "Ecos2" or "Acos2"'
      end select
      diff_cep=(phi_cep2-0.25d0)*2.d0-int((phi_cep2-0.25d0)*2.d0)
      if(ae_shape2=="Ecos2".and.abs(diff_cep)>=1.d-12)then
        stop "phi_cep2 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape2."
      end if
      if(rlaser_int_wcm2_2/=-1d0) &
        amplitude2=sqrt(rlaser_int_wcm2_2)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    
    !write information
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      if((fw%inc_dist1=='none').or.(fw%inc_dist2=='none')) then
        write(*,*) "Incident current source is placed at"
      else
        write(*,*) "Incident current sources are placed at"
      end if
      if(fw%inc_dist1/='none') then
        ii=1
        write(*,'(I8,A,3ES14.5,A)') ii,":",(fw%coo(fw%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        write(*,'(A,3ES14.5)') " ek_dir1:",ek_dir1
      end if
      if(fw%inc_dist2/='none') then
        ii=2
        write(*,'(I8,A,3ES14.5,A)') ii,":",(fw%coo(fw%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        write(*,'(A,3ES14.5)') " ek_dir2:",ek_dir2
      end if
      write(*,*) "**************************"
    end if
  end if
  
  !prepare linear response
  if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
    !check condition
    iflag=0
    if(iperiodic==3.and.trans_longi/='tr') iflag=1
    do ii=0,imedia_num
      if(fw%rep(ii)/=1.0d0.or.fw%rmu(ii)/=1.0d0.or.fw%sig(ii)/=0.0d0) iflag=1
      if(ii==0) then
        select case(type_media(ii))
        case('VACUUM','Vacuum','vacuum')
          continue
        case default
          iflag=1
        end select
      else
        select case(type_media(ii))
        case('DRUDE','Drude','drude','D','d')
          continue
        case default
          iflag=1
        end select
      end if
    end do
    if(iflag==1) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid input keywords:"
        write(*,*) "epsilon and rmu must be 1.0d0."
        write(*,*) "sigma must be 0.0d0."
        write(*,*) "type_media(i) must be drude, where i > 0."
        if(iperiodic==3) write(*,*) "trans_longi must be tr."
      end if
      stop
    end if
    
    !set initial current density
    if(fw%inum_d>0) then
      do ii=1,fw%inum_d
        do iz=fs%ng_sta(3),fs%ng_end(3)
        do iy=fs%ng_sta(2),fs%ng_end(2)
        do ix=fs%ng_sta(1),fs%ng_end(1)
          if(fw%idx_d(ix,iy,iz,ii)==1) then
            if(ae_shape1=='impulse') then
              fw%rjx_d(ix,iy,iz,ii)=fw%rjx_d(ix,iy,iz,ii) &
                                    -(omega_p_d(ii)**2.0d0)/(4.0d0*pi)*e_impulse*(epdir_re1(1)+epdir_im1(1))
            end if
            if(ae_shape2=='impulse') then
              fw%rjx_d(ix,iy,iz,ii)=fw%rjx_d(ix,iy,iz,ii) &
                                    -(omega_p_d(ii)**2.0d0)/(4.0d0*pi)*e_impulse*(epdir_re2(1)+epdir_im2(1))
            end if
          end if
          if(fw%idy_d(ix,iy,iz,ii)==1) then
            if(ae_shape1=='impulse') then
              fw%rjy_d(ix,iy,iz,ii)=fw%rjy_d(ix,iy,iz,ii) &
                                    -(omega_p_d(ii)**2.0d0)/(4.0d0*pi)*e_impulse*(epdir_re1(2)+epdir_im1(2))
            end if
            if(ae_shape2=='impulse') then
              fw%rjy_d(ix,iy,iz,ii)=fw%rjy_d(ix,iy,iz,ii) &
                                    -(omega_p_d(ii)**2.0d0)/(4.0d0*pi)*e_impulse*(epdir_re2(2)+epdir_im2(2))
            end if
          end if
          if(fw%idz_d(ix,iy,iz,ii)==1) then
            if(ae_shape1=='impulse') then
              fw%rjz_d(ix,iy,iz,ii)=fw%rjz_d(ix,iy,iz,ii) &
                                    -(omega_p_d(ii)**2.0d0)/(4.0d0*pi)*e_impulse*(epdir_re1(3)+epdir_im1(3))
            end if
            if(ae_shape2=='impulse') then
              fw%rjz_d(ix,iy,iz,ii)=fw%rjz_d(ix,iy,iz,ii) &
                                    -(omega_p_d(ii)**2.0d0)/(4.0d0*pi)*e_impulse*(epdir_re2(3)+epdir_im2(3))
            end if
          end if
        end do
        end do
        end do
      end do
    end if
    
    !initialize and allocate
    allocate(fw%time_lr(nt_em))
    fw%time_lr(:)=0.0d0
    fw%iter_lr=1
    allocate(fw%fr_lr(0:nenergy,3),fw%fi_lr(0:nenergy,3))
    fw%fr_lr(:,:)=0.0d0; fw%fi_lr(:,:)=0.0d0;
    allocate(fw%rjx_lr(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%rjy_lr(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             fw%rjz_lr(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                       fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                       fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd) )
    fw%rjx_lr(:,:,:)=0.0d0; fw%rjy_lr(:,:,:)=0.0d0; fw%rjz_lr(:,:,:)=0.0d0;
    if(iperiodic==0) then
      allocate(fw%px_lr(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                         fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                         fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
               fw%py_lr(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                         fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                         fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
               fw%pz_lr(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                         fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                         fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd) )
      fw%px_lr(:,:,:)=0.0d0; fw%py_lr(:,:,:)=0.0d0; fw%pz_lr(:,:,:)=0.0d0;
      allocate(fw%dip_lr(nt_em,3))
      fw%dip_lr(:,:)=0.0d0
    elseif(iperiodic==3) then
      allocate(fw%curr_lr(nt_em,3))
      fw%curr_lr(:,:)=0.0d0
    end if
  end if
  
  !write strat
  if(comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "**************************"
    write(*,*) "FDTD start"
    write(*,*) "**************************"
    write(*,*) "timestep"
    write(*,*) "-------------------------------------------------------"
  end if
  
contains
  
  !=========================================================================================
  != e and h allocation ====================================================================
  subroutine eh_allocate
    implicit none
    
    !e
    allocate(fw%ex_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_ex_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_ex_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ex_y(:,:,:)=0.0d0; fw%c1_ex_y(:,:,:)=0.0d0; fw%c2_ex_y(:,:,:)=0.0d0;
    allocate(fw%ex_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_ex_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_ex_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ex_z(:,:,:)=0.0d0; fw%c1_ex_z(:,:,:)=0.0d0; fw%c2_ex_z(:,:,:)=0.0d0;
    allocate(fw%ey_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_ey_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_ey_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ey_z(:,:,:)=0.0d0; fw%c1_ey_z(:,:,:)=0.0d0; fw%c2_ey_z(:,:,:)=0.0d0;
    allocate(fw%ey_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_ey_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_ey_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ey_x(:,:,:)=0.0d0; fw%c1_ey_x(:,:,:)=0.0d0; fw%c2_ey_x(:,:,:)=0.0d0;
    allocate(fw%ez_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_ez_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_ez_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ez_x(:,:,:)=0.0d0; fw%c1_ez_x(:,:,:)=0.0d0; fw%c2_ez_x(:,:,:)=0.0d0;
    allocate(fw%ez_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_ez_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_ez_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%ez_y(:,:,:)=0.0d0; fw%c1_ez_y(:,:,:)=0.0d0; fw%c2_ez_y(:,:,:)=0.0d0;
    
    !h
    allocate(fw%hx_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_hx_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_hx_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%hx_y(:,:,:)=0.0d0; fw%c1_hx_y(:,:,:)=0.0d0; fw%c2_hx_y(:,:,:)=0.0d0;
    allocate(fw%hx_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_hx_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_hx_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%hx_z(:,:,:)=0.0d0; fw%c1_hx_z(:,:,:)=0.0d0; fw%c2_hx_z(:,:,:)=0.0d0;
    allocate(fw%hy_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_hy_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_hy_z(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%hy_z(:,:,:)=0.0d0; fw%c1_hy_z(:,:,:)=0.0d0; fw%c2_hy_z(:,:,:)=0.0d0;
    allocate(fw%hy_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_hy_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_hy_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%hy_x(:,:,:)=0.0d0; fw%c1_hy_x(:,:,:)=0.0d0; fw%c2_hy_x(:,:,:)=0.0d0;
    allocate(fw%hz_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_hz_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_hz_x(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%hz_x(:,:,:)=0.0d0; fw%c1_hz_x(:,:,:)=0.0d0; fw%c2_hz_x(:,:,:)=0.0d0;
    allocate(fw%hz_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                     fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                     fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c1_hz_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    allocate(fw%c2_hz_y(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                        fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                        fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    fw%hz_y(:,:,:)=0.0d0; fw%c1_hz_y(:,:,:)=0.0d0; fw%c2_hz_y(:,:,:)=0.0d0;
    
  end subroutine eh_allocate
  
  !=========================================================================================
  != set fdtd coefficient ==================================================================
  subroutine eh_coeff
    implicit none
    real(8)  :: c1_e,c2_e_x,c2_e_y,c2_e_z,c1_h,c2_h_x,c2_h_y,c2_h_z,c2_j,&
                c1_e_mid,c2_e_x_mid,c2_e_y_mid,c2_e_z_mid,c2_j_mid
    
    !set constant parameter
    fw%rep(ii)=epsilon(ii); fw%rmu(ii)=rmu(ii); fw%sig(ii)=sigma(ii);
    
    !prepare coefficient
    c1_e  =(1.0d0-2.0d0*pi*fw%sig(ii)/fw%rep(ii)*fs%dt) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*fs%dt)
    c2_e_x=(fw%c_0/fw%rep(ii)*fs%dt) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*fs%dt)/fs%hgs(1)
    c2_e_y=(fw%c_0/fw%rep(ii)*fs%dt) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*fs%dt)/fs%hgs(2)
    c2_e_z=(fw%c_0/fw%rep(ii)*fs%dt) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*fs%dt)/fs%hgs(3)
    call comm_bcast(c1_e,  nproc_group_global)
    call comm_bcast(c2_e_x,nproc_group_global)
    call comm_bcast(c2_e_y,nproc_group_global)
    call comm_bcast(c2_e_z,nproc_group_global)
    c1_h=1.0d0
    c2_h_x=fw%c_0/fw%rmu(ii)*fs%dt/fs%hgs(1)
    c2_h_y=fw%c_0/fw%rmu(ii)*fs%dt/fs%hgs(2)
    c2_h_z=fw%c_0/fw%rmu(ii)*fs%dt/fs%hgs(3)
    call comm_bcast(c1_h,  nproc_group_global)
    call comm_bcast(c2_h_x,nproc_group_global)
    call comm_bcast(c2_h_y,nproc_group_global)
    call comm_bcast(c2_h_z,nproc_group_global)
    c2_j=(4.0d0*pi/fw%rep(ii)*fs%dt) &
         /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*fs%dt)
    call comm_bcast(c2_j,nproc_group_global)
    
    !check type_media
    select case(type_media(ii))
    case('PEC','Pec','pec')
      c1_e=0.0d0; c2_e_x=0.0d0; c2_e_y=0.0d0; c2_e_z=0.0d0;
    case('DRUDE','Drude','drude','D','d')
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        if(fs%imedia(ix,iy,iz)==ii) then
          if(fs%imedia(ix+1,iy,iz)==ii) then !x
            fw%idx_d(ix,iy,iz,icount_d)=1;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)<ii) then
            fw%idx_d(ix,iy,iz,icount_d)=1;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)>ii) then
            do ij=1,fw%inum_d
              if(fw%imedia_d(ij)==fs%imedia(ix+1,iy,iz)) then
                fw%idx_d(ix,iy,iz,ij)=1;
              end if
            end do
          end if
          if(fs%imedia(ix,iy+1,iz)==ii) then !y
            fw%idy_d(ix,iy,iz,icount_d)=1;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)<ii) then
            fw%idy_d(ix,iy,iz,icount_d)=1;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)>ii) then
            do ij=1,fw%inum_d
              if(fw%imedia_d(ij)==fs%imedia(ix,iy+1,iz)) then
                fw%idy_d(ix,iy,iz,ij)=1;
              end if
            end do
          end if
          if(fs%imedia(ix,iy,iz+1)==ii) then !z
            fw%idz_d(ix,iy,iz,icount_d)=1;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)<ii) then
            fw%idz_d(ix,iy,iz,icount_d)=1;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)>ii) then
            do ij=1,fw%inum_d
              if(fw%imedia_d(ij)==fs%imedia(ix,iy,iz+1)) then
                fw%idz_d(ix,iy,iz,ij)=1;
              end if
            end do
          end if
        end if
      end do
      end do
      end do
      fw%c1_j_d(icount_d)=(1.0d0-gamma_d(ii)*fs%dt/2.0d0)           / (1.0d0+gamma_d(ii)*fs%dt/2.0d0);
      fw%c2_j_d(icount_d)=((omega_p_d(ii)**2.0d0)*fs%dt/(4.0d0*pi)) / (1.0d0+gamma_d(ii)*fs%dt/2.0d0);
      icount_d=icount_d+1
    end select
    
    !set coefficient
    if(ii==0) then
      fw%c1_ex_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_e
      fw%c2_ex_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c2_e_y
      fw%c1_ex_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_e
      fw%c2_ex_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=-c2_e_z
          
      fw%c1_ey_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_e
      fw%c2_ey_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c2_e_z
      fw%c1_ey_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_e
      fw%c2_ey_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=-c2_e_x
        
      fw%c1_ez_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_e
      fw%c2_ez_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c2_e_x
      fw%c1_ez_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_e
      fw%c2_ez_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=-c2_e_y
        
      fw%c1_hx_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_h
      fw%c2_hx_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=-c2_h_y
      fw%c1_hx_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_h
      fw%c2_hx_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c2_h_z
        
      fw%c1_hy_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_h
      fw%c2_hy_z(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=-c2_h_z
      fw%c1_hy_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_h
      fw%c2_hy_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c2_h_x
      
      fw%c1_hz_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_h
      fw%c2_hz_x(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=-c2_h_x
      fw%c1_hz_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c1_h
      fw%c2_hz_y(fs%ng_sta(1):fs%ng_end(1),&
                 fs%ng_sta(2):fs%ng_end(2),&
                 fs%ng_sta(3):fs%ng_end(3))=c2_h_y
                  
      fw%c2_jx(fs%ng_sta(1):fs%ng_end(1),&
               fs%ng_sta(2):fs%ng_end(2),&
               fs%ng_sta(3):fs%ng_end(3))=-c2_j
      fw%c2_jy(fs%ng_sta(1):fs%ng_end(1),&
               fs%ng_sta(2):fs%ng_end(2),&
               fs%ng_sta(3):fs%ng_end(3))=-c2_j
      fw%c2_jz(fs%ng_sta(1):fs%ng_end(1),&
               fs%ng_sta(2):fs%ng_end(2),&
               fs%ng_sta(3):fs%ng_end(3))=-c2_j
    else
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        if(fs%imedia(ix,iy,iz)==ii) then
          !ex and jx
          if(fs%imedia(ix+1,iy,iz)==ii) then
            fw%c1_ex_y(ix,iy,iz)=c1_e; fw%c2_ex_y(ix,iy,iz)= c2_e_y;
            fw%c1_ex_z(ix,iy,iz)=c1_e; fw%c2_ex_z(ix,iy,iz)=-c2_e_z;
            fw%c2_jx(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)<ii) then
            fw%c1_ex_y(ix,iy,iz)=c1_e; fw%c2_ex_y(ix,iy,iz)= c2_e_y;
            fw%c1_ex_z(ix,iy,iz)=c1_e; fw%c2_ex_z(ix,iy,iz)=-c2_e_z;
            fw%c2_jx(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)>ii) then
            c1_e_mid  =(1.0d0-2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt)
            c2_e_y_mid=(fw%c_0/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /fs%hgs(2)
            c2_e_z_mid=(fw%c_0/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /fs%hgs(3)
            c2_j_mid  =(4.0d0*pi/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt)
            fw%c1_ex_y(ix,iy,iz)=c1_e_mid; fw%c2_ex_y(ix,iy,iz)= c2_e_y_mid;
            fw%c1_ex_z(ix,iy,iz)=c1_e_mid; fw%c2_ex_z(ix,iy,iz)=-c2_e_z_mid;
            fw%c2_jx(ix,iy,iz)=-c2_j_mid;
          end if
          
          !ey and jy
          if(fs%imedia(ix,iy+1,iz)==ii) then
            fw%c1_ey_z(ix,iy,iz)=c1_e; fw%c2_ey_z(ix,iy,iz)= c2_e_z;
            fw%c1_ey_x(ix,iy,iz)=c1_e; fw%c2_ey_x(ix,iy,iz)=-c2_e_x;
            fw%c2_jy(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)<ii) then
            fw%c1_ey_z(ix,iy,iz)=c1_e; fw%c2_ey_z(ix,iy,iz)= c2_e_z;
            fw%c1_ey_x(ix,iy,iz)=c1_e; fw%c2_ey_x(ix,iy,iz)=-c2_e_x;
            fw%c2_jy(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)>ii) then
            c1_e_mid  =(1.0d0-2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt)
            c2_e_z_mid=(fw%c_0/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt) &
                       /fs%hgs(3)
            c2_e_x_mid=(fw%c_0/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt) &
                       /fs%hgs(1)
            c2_j_mid  =(4.0d0*pi/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*fs%dt)
            fw%c1_ey_z(ix,iy,iz)=c1_e_mid; fw%c2_ey_z(ix,iy,iz)= c2_e_z_mid;
            fw%c1_ey_x(ix,iy,iz)=c1_e_mid; fw%c2_ey_x(ix,iy,iz)=-c2_e_x_mid;
            fw%c2_jy(ix,iy,iz)=-c2_j_mid;
          end if
          
          !ez and jz
          if(fs%imedia(ix,iy,iz+1)==ii) then
            fw%c1_ez_x(ix,iy,iz)=c1_e; fw%c2_ez_x(ix,iy,iz)= c2_e_x;
            fw%c1_ez_y(ix,iy,iz)=c1_e; fw%c2_ez_y(ix,iy,iz)=-c2_e_y;
            fw%c2_jz(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)<ii) then
            fw%c1_ez_x(ix,iy,iz)=c1_e; fw%c2_ez_x(ix,iy,iz)= c2_e_x;
            fw%c1_ez_y(ix,iy,iz)=c1_e; fw%c2_ez_y(ix,iy,iz)=-c2_e_y;
            fw%c2_jz(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)>ii) then
            c1_e_mid  =(1.0d0-2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*fs%dt)
            c2_e_x_mid=(fw%c_0/epsilon(fs%imedia(ix,iy,iz+1))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*fs%dt) &
                       /fs%hgs(1)
            c2_e_y_mid=(fw%c_0/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*fs%dt) &
                       /fs%hgs(2)
            c2_j_mid  =(4.0d0*pi/epsilon(fs%imedia(ix,iy,iz+1))*fs%dt) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*fs%dt)
            fw%c1_ez_x(ix,iy,iz)=c1_e_mid; fw%c2_ez_x(ix,iy,iz)= c2_e_x_mid;
            fw%c1_ez_y(ix,iy,iz)=c1_e_mid; fw%c2_ez_y(ix,iy,iz)=-c2_e_y_mid;
            fw%c2_jz(ix,iy,iz)=-c2_j_mid;
          end if
          
          !hx
          fw%c1_hx_y(ix,iy-1:iy,iz-1:iz)=c1_h; fw%c2_hx_y(ix,iy-1:iy,iz-1:iz)=-c2_h_y;
          fw%c1_hx_z(ix,iy-1:iy,iz-1:iz)=c1_h; fw%c2_hx_z(ix,iy-1:iy,iz-1:iz)= c2_h_z;
          
          !hy
          fw%c1_hy_z(ix-1:ix,iy,iz-1:iz)=c1_h; fw%c2_hy_z(ix-1:ix,iy,iz-1:iz)=-c2_h_z;
          fw%c1_hy_x(ix-1:ix,iy,iz-1:iz)=c1_h; fw%c2_hy_x(ix-1:ix,iy,iz-1:iz)= c2_h_x;
          
          !hz
          fw%c1_hz_x(ix-1:ix,iy-1:iy,iz)=c1_h; fw%c2_hz_x(ix-1:ix,iy-1:iy,iz)=-c2_h_x;
          fw%c1_hz_y(ix-1:ix,iy-1:iy,iz)=c1_h; fw%c2_hz_y(ix-1:ix,iy-1:iy,iz)= c2_h_y;
        end if
      end do
      end do
      end do
    end if
    
  end subroutine eh_coeff
  
  !=========================================================================================
  != apply smoothing =======================================================================
  subroutine eh_smoothing
    implicit none
    integer :: icomp
    
    if(fw%inum_d>0) then
      fw%wex_d(:,:,:,:)=dble(fw%idx_d(:,:,:,:))
      fw%wey_d(:,:,:,:)=dble(fw%idy_d(:,:,:,:))
      fw%wez_d(:,:,:,:)=dble(fw%idz_d(:,:,:,:))
      if(smooth_d=='y') then
        allocate(fs%imedia(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                           fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                           fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
        fs%imedia(:,:,:)=0
        allocate(fw%rmedia(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                           fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                           fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
        fw%rmedia(:,:,:)=0.0d0
        do ii=1,fw%inum_d
          do icomp=1,3
            if(icomp==1)     then
              fw%rmedia(:,:,:)=dble(fw%idx_d(:,:,:,ii))
            elseif(icomp==2) then
              fw%rmedia(:,:,:)=dble(fw%idy_d(:,:,:,ii))
            elseif(icomp==3) then
              fw%rmedia(:,:,:)=dble(fw%idz_d(:,:,:,ii))
            end if
            call eh_sendrecv(fs,fw,'r')
            fs%imedia(:,:,:)=int(fw%rmedia(:,:,:)+1d-3)
            do iz=fs%ng_sta(3),fs%ng_end(3)
            do iy=fs%ng_sta(2),fs%ng_end(2)
            do ix=fs%ng_sta(1),fs%ng_end(1)
              if(fs%imedia(ix,iy,iz)==1) then
                if(fs%imedia(ix+1,iy,iz)==0 .or. fs%imedia(ix-1,iy,iz)==0 .or. &
                   fs%imedia(ix,iy+1,iz)==0 .or. fs%imedia(ix,iy-1,iz)==0 .or. &
                   fs%imedia(ix,iy,iz+1)==0 .or. fs%imedia(ix,iy,iz-1)==0)then
                  if(icomp==1)     then
                    fw%wex_d(ix,iy,iz,ii)=weight_d
                  elseif(icomp==2) then
                    fw%wey_d(ix,iy,iz,ii)=weight_d
                  elseif(icomp==3) then
                    fw%wez_d(ix,iy,iz,ii)=weight_d
                  end if
                end if
              end if
            end do
            end do
            end do
          end do
        end do
        deallocate(fs%imedia); deallocate(fw%rmedia);
      end if
    end if
    
  end subroutine eh_smoothing
  
  !=========================================================================================
  != set pml ===============================================================================
  subroutine eh_set_pml(idir,c1_e1,c2_e1,c1_e2,c2_e2,c1_h1,c2_h1,c1_h2,c2_h2)
    implicit none
    integer,intent(in)  :: idir
    real(8),intent(out) :: c1_e1(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c2_e1(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c1_e2(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c2_e2(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c1_h1(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c2_h1(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c1_h2(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                           c2_h2(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                                 fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                                 fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd)
    integer :: ista,iend
    real(8) :: pml_del,s_max
    real(8) :: s_l(fw%ipml_l+1),sh_l(fw%ipml_l), &
               c1_pml(fw%ipml_l+1),c2_pml(fw%ipml_l+1),c1_pml_h(fw%ipml_l),c2_pml_h(fw%ipml_l)

    !set pml conductivity
    pml_del=fs%hgs(idir)
    s_max=-(fw%pml_m+1.0d0)*log(fw%pml_r)/(2.0d0*dble(fw%ipml_l)*pml_del) &
          *fw%c_0/(4.0d0*pi)*sqrt(fw%rep(0)/fw%rmu(0));
    do ii=1,(fw%ipml_l+1)
      s_l(ii)=s_max*(&
                    (dble(fw%ipml_l)*pml_del-(dble(ii)-1.0d0)*pml_del)/(dble(fw%ipml_l)*pml_del)&
                    )**fw%pml_m;
    end do
    do ii=1,fw%ipml_l
      sh_l(ii)=(fw%rmu(0)/fw%rep(0)) &
               *s_max*(&
                      (dble(fw%ipml_l)*pml_del-(dble(ii)-0.5d0)*pml_del)/(dble(fw%ipml_l)*pml_del)&
                      )**fw%pml_m;
    end do
    
    !set pml coefficient
    do ii=1,(fw%ipml_l+1)
      c1_pml(ii)=(1.0d0-2.0d0*pi*s_l(ii)/fw%rep(0)*fs%dt) &
                 /(1.0d0+2.0d0*pi*s_l(ii)/fw%rep(0)*fs%dt);
      c2_pml(ii)=(fw%c_0/fw%rep(0)*fs%dt) &
                 /(1.0d0+2.0d0*pi*s_l(ii)/fw%rep(0)*fs%dt)/pml_del
    end do
    call comm_bcast(c1_pml,nproc_group_global)
    call comm_bcast(c2_pml,nproc_group_global)
    do ii=1,fw%ipml_l
      c1_pml_h(ii)=(1.0d0-2.0d0*pi*sh_l(ii)/fw%rmu(0)*fs%dt) &
                   /(1.0d0+2.0d0*pi*sh_l(ii)/fw%rmu(0)*fs%dt);
      c2_pml_h(ii)=(fw%c_0/fw%rmu(0)*fs%dt) &
                   /(1.0d0+2.0d0*pi*sh_l(ii)/fw%rmu(0)*fs%dt)/pml_del
    end do
    call comm_bcast(c1_pml_h,nproc_group_global)
    call comm_bcast(c2_pml_h,nproc_group_global)
    
    !set pml(bottom)
    if((fs%i_bc(idir,1)==1).and.(fs%ng_sta(idir)<=(fs%lg_sta(idir)+fw%ipml_l))) then
      !e
      iend=fs%lg_sta(idir)+fw%ipml_l
      if(fs%ng_end(idir)<iend) then
        iend=fs%ng_end(idir)
      end if
      icount=1
      do ii=fs%ng_sta(idir),iend
        if(idir==1) then
          c1_e1(ii,:,:)= c1_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_e1(ii,:,:)=-c2_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c1_e2(ii,:,:)= c1_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_e2(ii,:,:)= c2_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
        elseif(idir==2) then
          c1_e1(:,ii,:)= c1_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_e1(:,ii,:)=-c2_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c1_e2(:,ii,:)= c1_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_e2(:,ii,:)= c2_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
        elseif(idir==3) then
          c1_e1(:,:,ii)= c1_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_e1(:,:,ii)=-c2_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c1_e2(:,:,ii)= c1_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_e2(:,:,ii)= c2_pml(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
        end if
        icount=icount+1
      end do
      
      !h
      if(iend==(fs%lg_sta(idir)+fw%ipml_l)) then
        iend=iend-1
      end if
      icount=1
      do ii=fs%ng_sta(idir),iend
        if(idir==1) then
          c1_h1(ii,:,:)= c1_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_h1(ii,:,:)= c2_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c1_h2(ii,:,:)= c1_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_h2(ii,:,:)=-c2_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
        elseif(idir==2) then
          c1_h1(:,ii,:)= c1_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_h1(:,ii,:)= c2_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c1_h2(:,ii,:)= c1_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_h2(:,ii,:)=-c2_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
        elseif(idir==3) then
          c1_h1(:,:,ii)= c1_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_h1(:,:,ii)= c2_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c1_h2(:,:,ii)= c1_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
          c2_h2(:,:,ii)=-c2_pml_h(fs%ng_sta(idir)-fs%lg_sta(idir)+icount)
        end if
        icount=icount+1
      end do
    end if
    
    !set pml(top)
    if((fs%i_bc(idir,2)==1).and.(fs%ng_end(idir)>=(fs%lg_end(idir)-fw%ipml_l))) then
      !e
      ista=fs%lg_end(idir)-fw%ipml_l
      if(fs%ng_sta(idir)>ista) then
        ista=fs%ng_sta(idir)
      end if
      icount=1
      do ii=ista,fs%ng_end(idir)
        if(idir==1) then
          c1_e1(ii,:,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_e1(ii,:,:)=-c2_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c1_e2(ii,:,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_e2(ii,:,:)= c2_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==2) then
          c1_e1(:,ii,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_e1(:,ii,:)=-c2_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c1_e2(:,ii,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_e2(:,ii,:)= c2_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==3) then
          c1_e1(:,:,ii)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_e1(:,:,ii)=-c2_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c1_e2(:,:,ii)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_e2(:,:,ii)= c2_pml((fw%ipml_l+1)-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
        end if
        icount=icount+1
      end do
      
      !h
      if(fs%ng_end(idir)==fs%lg_end(idir)) then
        iend=fs%ng_end(idir)-1
      else
        iend=fs%ng_end(idir)
      end if
      icount=1
      do ii=ista,iend
        if(idir==1) then
          c1_h1(ii,:,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_h1(ii,:,:)= c2_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c1_h2(ii,:,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_h2(ii,:,:)=-c2_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==2) then
          c1_h1(:,ii,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_h1(:,ii,:)= c2_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c1_h2(:,ii,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_h2(:,ii,:)=-c2_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==3) then
          c1_h1(:,:,ii)= c1_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_h1(:,:,ii)= c2_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c1_h2(:,:,ii)= c1_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
          c2_h2(:,:,ii)=-c2_pml_h(fw%ipml_l-(ista-(fs%lg_end(idir)-fw%ipml_l)+(icount-1)))
        end if
        icount=icount+1
      end do
    end if
    
  end subroutine eh_set_pml
  
end subroutine eh_init

!=========================================================================================
!= input fdtd shape data =================================================================
subroutine eh_input_shape(ifn,ng_sta,ng_end,lg_sta,lg_end,Nd,imat,format)
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use inputoutput,          only: shape_file
  implicit none
  integer,intent(in)      :: ifn,Nd
  integer,intent(in)      :: ng_sta(3),ng_end(3),lg_sta(3),lg_end(3)
  integer,intent(out)     :: imat(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                  ng_sta(2)-Nd:ng_end(2)+Nd,&
                                  ng_sta(3)-Nd:ng_end(3)+Nd)
  character(2),intent(in) :: format
  real(8),allocatable     :: rtmp1d(:)
  integer                 :: inum(3),inum_check(3)
  integer                 :: ii,ij,ix,iy,iz,iflag_x,iflag_y,iflag_z
  
  !open file
  open(ifn,file=trim(shape_file), status='old')
  
  if(trim(format)=='cu') then
    !check grid information
    inum(:)=lg_end(:)-lg_sta(:)+1
    read(ifn,*); read (ifn,*); read (ifn,*); !skip
    allocate(rtmp1d(4))
    read (ifn,*) rtmp1d; inum_check(1)=int(rtmp1d(1)+1d-3);
    read (ifn,*) rtmp1d; inum_check(2)=int(rtmp1d(1)+1d-3);
    read (ifn,*) rtmp1d; inum_check(3)=int(rtmp1d(1)+1d-3);
    deallocate(rtmp1d)
    do ii=1,3
      if(inum(ii)/=inum_check(ii)) then
        if(comm_is_root(nproc_id_global)) write(*,*) "al_em or dl_em does not mutch shape file."
        stop
      end if
    end do
    read (ifn,*); !skip
    
    !input shape(general case)
    allocate(rtmp1d(6))
    ix=lg_sta(1); iy=lg_sta(2); iz=lg_sta(3);
    do ii=1,int(inum(1)*inum(2)*inum(3)/6)
      read (ifn,*) rtmp1d
      do ij=1,6
        !check flag and write imat
        iflag_x=0; iflag_y=0; iflag_z=0;
        if(ix>=ng_sta(1) .and. ix<=ng_end(1)) iflag_x=1
        if(iy>=ng_sta(2) .and. iy<=ng_end(2)) iflag_y=1
        if(iz>=ng_sta(3) .and. iz<=ng_end(3)) iflag_z=1
        if(iflag_x==1 .and. iflag_y==1 .and. iflag_z==1) then
          imat(ix,iy,iz)=int(rtmp1d(ij)+1d-3)
        end if        
        
        !update iz, iy, ix 
        iz=iz+1                                            !iz
        if(iz>lg_end(3))                      iz=lg_sta(3) !iz
        if(iz==lg_sta(3))                     iy=iy+1      !iy
        if(iy>lg_end(2))                      iy=lg_sta(2) !iy
        if(iz==lg_sta(3) .and. iy==lg_sta(2)) ix=ix+1      !ix
      end do
    end do
    deallocate(rtmp1d)
    
    !input shape(special case)
    if(mod(inum(1)*inum(2)*inum(3),6)>0) then
      allocate(rtmp1d(mod(inum(1)*inum(2)*inum(3),6)))
      read (ifn,*) rtmp1d
      do ij=1,mod(inum(1)*inum(2)*inum(3),6)
        !check flag and write imat
        iflag_x=0; iflag_y=0; iflag_z=0;
        if(ix>=ng_sta(1) .and. ix<=ng_end(1)) iflag_x=1
        if(iy>=ng_sta(2) .and. iy<=ng_end(2)) iflag_y=1
        if(iz>=ng_sta(3) .and. iz<=ng_end(3)) iflag_z=1
        if(iflag_x==1 .and. iflag_y==1 .and. iflag_z==1) then
          imat(ix,iy,iz)=int(rtmp1d(ij)+1d-3)
        end if        
        
        !update iz, iy, ix 
        iz=iz+1                                            !iz
        if(iz>lg_end(3))                      iz=lg_sta(3) !iz
        if(iz==lg_sta(3))                     iy=iy+1      !iy
        if(iy>lg_end(2))                      iy=lg_sta(2) !iy
        if(iz==lg_sta(3) .and. iy==lg_sta(2)) ix=ix+1      !ix
      end do      
      deallocate(rtmp1d)
    end if
  elseif(trim(format)=='mp') then
  end if
  
  !close file
  close(ifn)

end subroutine eh_input_shape

!=========================================================================================
!= find point and set corresponding processor element ====================================
subroutine eh_find_point(rloc,id,ipo,ili,ipl,ista,iend,icoo_sta,icoo_end,coo)
  use salmon_parallel,      only: nproc_id_global,nproc_size_global,nproc_group_global
  use salmon_communication, only: comm_summation
  implicit none
  real(8),intent(in)    :: rloc(3)
  integer,intent(inout) :: id(3)
  integer,intent(out)   :: ipo
  integer,intent(out)   :: ili(3),ipl(3)
  integer,intent(in)    :: ista(3),iend(3)
  integer,intent(in)    :: icoo_sta,icoo_end
  real(8),intent(in)    :: coo(icoo_sta:icoo_end,3)
  integer               :: ii,ix,iy,iz,ipe_tmp,i1,i1s,i2,i2s
  integer               :: id_tmp(3),id_tmp2(3)
  real(8)               :: err(0:nproc_size_global-1),err2(0:nproc_size_global-1)
  real(8)               :: err_tmp
  
  !set initial value
  err(:)              =0.0d0
  err(nproc_id_global)=1.0d9
  err2(:)             =0.0d0
  id_tmp(:)           =0
  id_tmp2(:)          =0
  
  !find observation point in each PE
  do iz=ista(3),iend(3)
  do iy=ista(2),iend(2)
  do ix=ista(1),iend(1)
    err_tmp=sqrt( (coo(ix,1)-rloc(1))**2.0d0 &
                 +(coo(iy,2)-rloc(2))**2.0d0 &
                 +(coo(iz,3)-rloc(3))**2.0d0 )
    if(err(nproc_id_global)>err_tmp) then
      err(nproc_id_global)=err_tmp
      id_tmp(1)=ix; id_tmp(2)=iy; id_tmp(3)=iz;
    end if
  end do
  end do
  end do
  call comm_summation(err,err2,nproc_size_global,nproc_group_global)
  
  !determine and share id + determine pe including the point
  ipe_tmp=-1; err_tmp=1.0d9;
  do ii=0,nproc_size_global-1
    if(err_tmp>err2(ii)) then
      err_tmp=err2(ii)
      ipe_tmp=ii
    end if
  end do
  if(nproc_id_global==ipe_tmp) then
    ipo=1;
  else
    ipo=0; id_tmp(:)=0;
  end if
  call comm_summation(id_tmp,id_tmp2,3,nproc_group_global)
  id(:)=id_tmp2(:)
  
  !determine pe including the line
  do ii=1,3
    if(ii==1) then     !x-line(searching at yz-plane)
      i1s=3; i2s=2;
    elseif(ii==2) then !x-line(searching at xz-plane)
      i1s=3; i2s=1;
    elseif(ii==3) then !z-line(searching at xy-plane)
      i1s=2; i2s=1;
    end if
    do i2=ista(i2s),iend(i2s)
    do i1=ista(i1s),iend(i1s)
      if( (i1==id(i1s)).and.(i2==id(i2s)) ) ili(ii)=1
    end do
    end do
  end do
  
  !determine pe including the plane
  do ii=1,3
    if(ii==1) then     !xy-plane(searching at z-line)
      i1s=3;
    elseif(ii==2) then !yz-plane(searching at x-line)
      i1s=1;
    elseif(ii==3) then !xz-plane(searching at y-line)
      i1s=2;
    end if
    do i1=ista(i1s),iend(i1s)
      if(i1==id(i1s)) ipl(ii)=1
    end do
  end do
  
end subroutine eh_find_point

!=========================================================================================
!= set coordinate ========================================================================
subroutine eh_set_coo(iperi,Nd,ioe,ista,iend,hgs,coo)
  implicit none
  integer,intent(in)  :: iperi,Nd
  integer,intent(in)  :: ioe(3),ista(3),iend(3)
  real(8),intent(in)  :: hgs(3)
  real(8),intent(out) :: coo(minval(ista(:))-Nd:maxval(iend(:))+Nd,3)
  integer :: ii,ij
  
  do ii=1,3
    select case(iperi)
    case(0)
      select case(ioe(ii))
      case(1)
        do ij=ista(ii)-Nd,iend(ii)+Nd
          coo(ij,ii)=dble(ij)*hgs(ii)
        end do
      case(2)
        do ij=ista(ii)-Nd,iend(ii)+Nd
          coo(ij,ii)=(dble(ij)-0.5d0)*hgs(ii)
        end do
      end select
    case(3)
      do ij=ista(ii)-Nd,iend(ii)+Nd
        coo(ij,ii)=dble(ij-1)*hgs(ii)
      end do
    end select
  end do
  
end subroutine eh_set_coo

!=========================================================================================
!= prepare GCEED =========================================================================
!= (This routine is temporary) ===========================================================
!= (With unifying ARTED and GCEED, this routine will be removed) =========================
subroutine eh_prep_GCEED(fs,fw)
  use inputoutput,       only: nproc_domain,nproc_domain_s,num_kgrid,nproc_k,nproc_ob,isequential,iperiodic
  use salmon_parallel,   only: nproc_id_orbitalgrid,nproc_id_global,nproc_size_global,nproc_group_global
  use set_numcpu,        only: set_numcpu_gs
  use scf_data,          only: nproc_Mxin,nproc_Mxin_s,nproc_Mxin_mul,nproc_Mxin_mul_s_dm,nproc_Mxin_s_dm,&
                               k_sta,k_end,k_num,num_kpoints_3d,num_kpoints_rd,&
                               rLsize,Harray,Hgs,Hvol,imesh_oddeven,&
                               lg_sta,lg_end,lg_num, &
                               mg_sta,mg_end,mg_num, &
                               ng_sta,ng_end,ng_num,&
                               ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,Nd, &
                               ista_Mxin,iend_Mxin,inum_Mxin,&
                               ista_Mxin_s,iend_Mxin_s,inum_Mxin_s
  use new_world_sub,     only: make_new_world
  use init_sendrecv_sub, only: init_updown,iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
  use sendrecv_grid,     only: init_sendrecv_grid
  use structures,        only: s_fdtd_system
  use salmon_maxwell,    only: s_fdtd_work
  implicit none
  type(s_fdtd_system)   :: fs
  type(s_fdtd_work)     :: fw
  integer               :: neig_ng_eh(1:3,1:2)
  integer               :: id_tmp,ii
  
  !set mpi condition
  num_kpoints_3d(1:3)=num_kgrid(1:3)
  num_kpoints_rd=num_kpoints_3d(1)*num_kpoints_3d(2)*num_kpoints_3d(3)
  nproc_Mxin=nproc_domain
  nproc_Mxin_s=nproc_domain_s
  call set_numcpu_gs(nproc_mxin,nproc_mxin_s,nproc_mxin_s_dm)
  nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
  nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  call make_new_world
  call setk(k_sta,k_end,k_num,num_kpoints_rd,nproc_k,nproc_id_orbitalgrid)
  
  !set grid and odd or even grid paterns
  rLsize(:,1)=fs%rlsize(:); Harray(:,1)=fs%hgs(:);
  Hgs(:)=Harray(:,1); Hvol=Hgs(1)*Hgs(2)*Hgs(3);
  call set_imesh_oddeven(1)
  call setlg(fs%lg,lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
             Hgs,Nd,rLsize(:,1),imesh_oddeven,iperiodic,1)
  allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1), &
           inum_Mxin(3,0:nproc_size_global-1))
  call setmg(fs%mg,mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
             lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_Mxin,nproc_k,nproc_ob,isequential,1)
  allocate(ista_Mxin_s(3,0:nproc_size_global-1),iend_Mxin_s(3,0:nproc_size_global-1))
  allocate(inum_Mxin_s(3,0:nproc_size_global-1))
  call setng(fs%ng,ng_sta,ng_end,ng_num,ista_Mxin_s,iend_Mxin_s,inum_Mxin_s, &
             nproc_size_global,nproc_id_global,nproc_Mxin,nproc_Mxin_s_dm,ista_Mxin,iend_Mxin,isequential,1)
  fs%lg_sta(:)=fs%lg%is(:); fs%lg_end(:)=fs%lg%ie(:);
  fs%ng_sta(:)=fs%ng%is(:); fs%ng_end(:)=fs%ng%ie(:);
  fw%ioddeven(:)=imesh_oddeven(:);
  
  !set sendrecv environment
  call init_updown
  if (iperiodic==0) then
    id_tmp=2;
  elseif (iperiodic==3) then
    !This process is temporal. 
    !With bug-fixing init_updown for iperiodic=3 and ob=1, this process will be removed.
    id_tmp=1;
  end if
  neig_ng_eh(1,1)=iup_array(id_tmp); neig_ng_eh(1,2)=idw_array(id_tmp);
  neig_ng_eh(2,1)=jup_array(id_tmp); neig_ng_eh(2,2)=jdw_array(id_tmp);
  neig_ng_eh(3,1)=kup_array(id_tmp); neig_ng_eh(3,2)=kdw_array(id_tmp);
  !This process about ng is temporal. 
  !With modifying set_ng to be applied to arbitrary Nd, this process will be removed.
  fs%ng%is_overlap(1:3)=fs%ng_sta(1:3)-fw%Nd
  fs%ng%ie_overlap(1:3)=fs%ng_end(1:3)+fw%Nd
  fs%ng%is_array(1:3)  =fs%ng_sta(1:3)-fw%Nd
  fs%ng%ie_array(1:3)  =fs%ng_end(1:3)+fw%Nd
  if(allocated(fs%ng%idx)) deallocate(fs%ng%idx)
  if(allocated(fs%ng%idy)) deallocate(fs%ng%idy)
  if(allocated(fs%ng%idz)) deallocate(fs%ng%idz)
  allocate(fs%ng%idx(fs%ng%is_overlap(1):fs%ng%ie_overlap(1)), &
           fs%ng%idy(fs%ng%is_overlap(2):fs%ng%ie_overlap(2)), &
           fs%ng%idz(fs%ng%is_overlap(3):fs%ng%ie_overlap(3)))
  do ii=fs%ng%is_overlap(1),fs%ng%ie_overlap(1)
    fs%ng%idx(ii)=ii
  end do
  do ii=fs%ng%is_overlap(2),fs%ng%ie_overlap(2)
    fs%ng%idy(ii)=ii
  end do
  do ii=fs%ng%is_overlap(3),fs%ng%ie_overlap(3)
    fs%ng%idz(ii)=ii
  end do
  fs%ng%nd=fw%Nd
  call init_sendrecv_grid(fs%srg_ng,fs%ng,1,nproc_group_global,nproc_id_global,neig_ng_eh)
  
end subroutine eh_prep_GCEED
