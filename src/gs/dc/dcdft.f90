!
!  Copyright 2019-2024 SALMON developers
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
module dcdft
  implicit none
contains

  subroutine init_dcdft(dc,mixing)
    use structures
    use parallelization, only: nproc_group_global
    use salmon_global, only: num_fragment, base_directory, nproc_k, nproc_ob, nproc_rgrid, &
    & nproc_rgrid_tot, nelec, nstate, yn_dc
    implicit none
    type(s_dcdft) ,intent(inout) :: dc
    type(s_mixing),intent(inout) :: mixing
    !
    integer :: i
    integer :: nproc_ob_tmp, nproc_rgrid_tmp(3)
    type(s_stencil) :: stencil_dummy
    type(s_sendrecv_grid):: srg_dummy
    type(s_ofile) :: ofile_dummy
    
    dc%icomm_tot = nproc_group_global
    dc%n_frag = num_fragment(1)*num_fragment(2)*num_fragment(3)
    dc%elec_num_tot = dble(nelec)
    dc%base_directory = trim(base_directory)
    
    dc%nstate_frag = nstate ! nstate for the fragment !!!!!! future work: new input variable
    
    if(nproc_k/=1) then
      stop "DC method (yn_dc=y): nproc_k must be 1 for both the total system and fragments."
    end if
    nproc_ob_tmp = nproc_ob
    nproc_rgrid_tmp = nproc_rgrid
    
  ! total system
    nproc_ob = 1 ! override
    nproc_rgrid = nproc_rgrid_tot ! override
    nstate = nelec ! override !!!!!! future work: remove
    call init_total
    
  ! fragment
    nproc_ob = nproc_ob_tmp ! override
    nproc_rgrid = nproc_rgrid_tmp ! override
    nstate = dc%nstate_frag ! override !!!!!! future work: new input variable
    call init_comm_frag
    call init_fragment
    
  contains
  
    subroutine init_total
      use initialization_sub, only: init_dft
      use sendrecv_grid, only: dealloc_cache
      use mixing_sub, only: init_mixing
      use salmon_pp, only: read_pslfile
      use prep_pp_sub, only: init_ps
      implicit none
      type(s_pp_info) :: pp_tmp
      type(s_pp_grid) :: ppg_tmp
      
      call init_dft(dc%icomm_tot,dc%info_tot,dc%lg_tot,dc%mg_tot,dc%system_tot, &
      & stencil_dummy,dc%fg_tot,dc%poisson_tot,srg_dummy,dc%srg_scalar_tot,ofile_dummy)
      deallocate(dc%system_tot%rocc)
      call dealloc_cache(srg_dummy)
      
      call allocate_scalar(dc%mg_tot,dc%rho_tot)
      call allocate_scalar(dc%mg_tot,dc%vh_tot)
      call allocate_scalar(dc%mg_tot,dc%vpsl_tot)
      allocate(dc%rho_tot_s(dc%system_tot%nspin),dc%vloc_tot(dc%system_tot%nspin),dc%vxc_tot(dc%system_tot%nspin))
      do i=1,dc%system_tot%nspin
        call allocate_scalar(dc%mg_tot,dc%rho_tot_s(i))
        call allocate_scalar(dc%mg_tot,dc%vloc_tot(i))
        call allocate_scalar(dc%mg_tot,dc%vxc_tot(i))
      end do
      
    ! mixing
      mixing%num_rho_stock = 21
      call init_mixing(dc%system_tot%nspin,dc%mg_tot,mixing)
      
    ! Vpsl
      call read_pslfile(dc%system_tot,pp_tmp)
      call init_ps(dc%lg_tot,dc%mg_tot,dc%system_tot,dc%info_tot,dc%fg_tot,dc%poisson_tot, &
      & pp_tmp,ppg_tmp,dc%vpsl_tot)
    
    end subroutine init_total
  
    subroutine init_comm_frag
      use parallelization
      use mpi !!!!! future work --> wrapper
      implicit none
      integer :: nproc,myrank
      integer :: comm_F,nproc_F,myrank_F
      integer :: npg,i,j,k,m,ierr
      
      !dc%icomm_tot = nproc_group_global ! (already)
      myrank = nproc_id_global
      nproc = nproc_size_global
      
      ! set dc%i_frag
      npg = nproc / dc%n_frag
      m = mod(nproc,dc%n_frag) ! nproc = npg*dc%n_frag + m
      k=0
      do j=0,dc%n_frag-1
      do i=0,npg-1
        if(j*npg+i==myrank) then
          dc%i_frag=j
          k=1
          exit
          exit
        end if
      end do
      end do
      if(k==0) dc%i_frag = myrank-npg*dc%n_frag
      dc%i_frag = dc%i_frag + 1 ! = 1:dc%n_frag
      
      ! split communicator
      call Mpi_Comm_split(dc%icomm_tot,dc%i_frag,myrank,comm_F,ierr) ! dc%i_frag : color, myrank : key
      call MPi_Comm_SIZE(comm_F,nproc_F,ierr)
      call mpi_comm_rank(comm_F,myrank_F,ierr)
      
      ! Override global variables
      nproc_group_global = comm_F
      nproc_id_global = myrank_F
      nproc_size_global = nproc_F
      
write(*,*) "test_dcdft 1:",myrank_F,nproc_F,dc%i_frag,myrank,nproc  !!!!!!!!! test_dcdft
    
    end subroutine init_comm_frag
    
    subroutine init_fragment
      use parallelization !!!!!!!!! test_dcdft
      use salmon_global, only: length_buffer, kion, rion, natom, num_rgrid, al
      implicit none
      integer :: i,j,k,n
      integer :: iatom,iatom_frag,nxyz_buffer(3)
      integer :: kion_frag(natom,dc%n_frag),natom_frag(dc%n_frag)
      real(8) :: dr
      real(8) :: r1(3),r2(3)
      real(8) :: ldomain(3)
      real(8) :: rion_frag(3,natom,dc%n_frag)
    
    ! length of domain
      ldomain(1:3) = al(1:3) / dble(num_fragment(1:3))
      
      do n=1,3 ! x,y,z
      ! rion --> rion = [0:al] (total system)
        do i=1,natom
          do while (rion(n,i) < 0)
            rion(n,i) = rion(n,i) + al(n)
          end do
          do while (rion(n,i) > al(n))
            rion(n,i) = rion(n,i) - al(n)
          end do
          if(rion(n,i) < 0 .or. rion(n,i) > al(n)) stop "DC method (yn_dc=y): rion"
        end do
      ! dc%nxyz_domain: # of grid points for each domain
        if(mod(num_rgrid(n),num_fragment(n))==0) then
          dc%nxyz_domain(n) = num_rgrid(n) / num_fragment(n)
          dr = al(n)/dble(num_rgrid(n))
          dr = length_buffer(n)/dr
          if(dr-dble(int(dr)) /= 0d0) stop "DC method (yn_dc=y): grid mismatch of length_buffer"
          nxyz_buffer(n) = int(dr)
        else
          stop "DC method (yn_dc=y): mod(num_rgrid,num_fragment) /= 0"
        end if
      end do

    ! position of each fragment
      allocate(dc%ixyz_frag(3,dc%n_frag),dc%rxyz_frag(3,dc%n_frag))
      n = 1
      do i=1,num_fragment(1)
      do j=1,num_fragment(2)
      do k=1,num_fragment(3)
      ! ix_total = ix_fragment + dc%ixyz_frag(1,dc%i_frag), ix_fragment=[1:dc%nxyz_domain(1)], etc.
        dc%ixyz_frag(1,n) = (i-1)*dc%nxyz_domain(1)
        dc%ixyz_frag(2,n) = (j-1)*dc%nxyz_domain(2)
        dc%ixyz_frag(3,n) = (k-1)*dc%nxyz_domain(3)
        dc%rxyz_frag(1,n) = dble(i-1)*ldomain(1)
        dc%rxyz_frag(2,n) = dble(j-1)*ldomain(2)
        dc%rxyz_frag(3,n) = dble(k-1)*ldomain(3)
        n = n + 1
      end do
      end do
      end do
      
    ! variables for each fragment
      n = 1
      do i=1,num_fragment(1)
      do j=1,num_fragment(2)
      do k=1,num_fragment(3)
        r1 = dc%rxyz_frag(:,n) - length_buffer
        r2 = r1 + ldomain + length_buffer
      ! atom count
        iatom_frag = 0
        do iatom=1,natom
          if(rion(1,iatom)>=r1(1) .and. rion(1,iatom)<r2(1) .and. &
          &  rion(2,iatom)>=r1(2) .and. rion(2,iatom)<r2(2) .and. &
          &  rion(3,iatom)>=r1(3) .and. rion(3,iatom)<r2(3) ) then
            iatom_frag = iatom_frag + 1
            rion_frag(1:3,iatom_frag,n) = rion(1:3,iatom) - dc%rxyz_frag(1:3,n)
            kion_frag(iatom_frag,n) = kion(iatom)
          end if
        end do
        natom_frag(n) = iatom_frag
        n = n + 1
      end do
      end do
      end do
    
    ! set variables for own fragment
    
    ! nelec (total system) --> nelec (fragment)
      nelec = nelec * natom_frag(dc%i_frag) / natom !!!!!!!!! test_dcdft !!!!!!!!
    
    ! al, num_rgrid (total system) --> al, num_rgrid (fragment)
      al = ldomain + 2d0*length_buffer
      num_rgrid = dc%nxyz_domain + 2*nxyz_buffer
      
    ! natom, rion, kion (total system) --> natom, rion, kion (fragment)
      natom = natom_frag(dc%i_frag)
      deallocate(rion,kion)
      allocate(rion(3,natom),kion(natom))
      rion(1:3,1:natom) = rion_frag(1:3,1:natom,dc%i_frag)
      kion(1:natom) = kion_frag(1:natom,dc%i_frag)
      
    ! base_directory =  !!!!!!!!!! future work
      
    ! dc%jxyz_tot: r-grid (fragment) --> r-grid (total)
      allocate(dc%jxyz_tot(maxval(num_rgrid),3))
      do n=1,3 ! x,y,z
        do i=1,num_rgrid(n) ! r-grid (fragment)
          if(i <= dc%nxyz_domain(n) + nxyz_buffer(n)) then
            j = i + dc%ixyz_frag(n,dc%i_frag)
          else
            j = ( i - num_rgrid(n) ) + dc%ixyz_frag(n,dc%i_frag) ! minus region
          end if
          j = mod(j+dc%lg_tot%num(n)-1,dc%lg_tot%num(n))+1
          dc%jxyz_tot(i,n) = j ! r-grid (total)
        end do
      end do
      
if(nproc_id_global==0) write(*,*) "test_dcdft 2: i_frag, natom, nelec",dc%i_frag, natom, nelec  !!!!!!!!! test_dcdft
    
    end subroutine init_fragment
    
  end subroutine init_dcdft
  
  ! rho_s (fragment) --> dc%rho_tot_s (total system) !!!!!! future work: occupancy, spsi
  subroutine calc_rho_total_dcdft(nspin,lg,mg,info,rho_s,dc)
    use structures
    use communication, only: comm_summation
    implicit none
    integer,              intent(in) :: nspin
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_parallel_info),intent(in) :: info
    type(s_scalar),       intent(in) :: rho_s(nspin)
    type(s_dcdft)                    :: dc
    !
    integer :: ix,iy,iz,ispin,ix_tot,iy_tot,iz_tot
    real(8),dimension(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:nspin) :: frg_tmp,frg
    real(8),dimension(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1:nspin) :: tot_tmp,tot
    
    ! rho_s (fragment)
    frg_tmp = 0d0
    do ispin=1,nspin
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      frg_tmp(ix,iy,iz,ispin) = rho_s(ispin)%f(ix,iy,iz)
    end do
    end do
    end do
    end do
    call comm_summation(frg_tmp,frg,lg%num(1)*lg%num(2)*lg%num(3)*nspin,info%icomm_r)
    
    ! rho_s (total)
    tot_tmp = 0d0
    if(info%id_rko==0) then
      do ispin=1,nspin
      do iz=1,dc%nxyz_domain(3); iz_tot = dc%ixyz_frag(3,dc%i_frag) + iz
      do iy=1,dc%nxyz_domain(2); iy_tot = dc%ixyz_frag(2,dc%i_frag) + iy
      do ix=1,dc%nxyz_domain(1); ix_tot = dc%ixyz_frag(1,dc%i_frag) + ix
        tot_tmp(ix_tot,iy_tot,iz_tot,ispin) = frg(ix,iy,iz,ispin)
      end do
      end do
      end do
      end do
    end if
    call comm_summation(tot_tmp,tot,dc%lg_tot%num(1)*dc%lg_tot%num(2)*dc%lg_tot%num(3)*nspin,dc%icomm_tot)
    do ispin=1,nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      dc%rho_tot_s(ispin)%f(ix,iy,iz) = tot(ix,iy,iz,ispin)
    end do
    end do
    end do
    end do
    
  end subroutine calc_rho_total_dcdft
  
  ! dc%vloc_tot (total system) --> v_local (fragment)
  subroutine calc_vlocal_fragment_dcdft(nspin,lg,mg,info,vloc,dc)
    use structures
    use communication, only: comm_summation
    implicit none
    integer,              intent(in) :: nspin
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_parallel_info),intent(in) :: info
    type(s_scalar)                   :: vloc(nspin)
    type(s_dcdft)                    :: dc
    !
    integer :: ix,iy,iz,ispin,ix_tot,iy_tot,iz_tot
    real(8),dimension(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),1:nspin) :: tot_tmp,tot
    
    ! vloc (total)
    tot_tmp = 0d0
    do ispin=1,nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      tot_tmp(ix,iy,iz,ispin) = dc%vloc_tot(ispin)%f(ix,iy,iz)
    end do
    end do
    end do
    end do
    call comm_summation(tot_tmp,tot,dc%lg_tot%num(1)*dc%lg_tot%num(2)*dc%lg_tot%num(3)*nspin,dc%icomm_tot)
    
    ! vloc (fragment)
    do ispin=1,nspin
    do iz=mg%is(3),mg%ie(3) ; iz_tot = dc%jxyz_tot(iz,3)
    do iy=mg%is(2),mg%ie(2) ; iy_tot = dc%jxyz_tot(iy,2)
    do ix=mg%is(1),mg%ie(1) ; ix_tot = dc%jxyz_tot(ix,1)
      vloc(ispin)%f(ix,iy,iz) = tot(ix_tot,iy_tot,iz_tot,ispin)
    end do
    end do
    end do
    end do
    
  end subroutine calc_vlocal_fragment_dcdft
  
  subroutine test_density(dc) !!!!!!! test_dcdft
    use structures
    use communication, only: comm_summation
    use salmon_global, only: natom, kion, rion
    use writefield, only: write_dns
    implicit none
    type(s_dcdft) :: dc
    !
    natom = dc%system_tot%nion
    deallocate(kion,rion)
    allocate(kion(natom),rion(3,natom))
    kion = dc%system_tot%kion
    rion = dc%system_tot%rion
    call write_dns(dc%lg_tot,dc%mg_tot,dc%system_tot,dc%info_tot,dc%rho_tot_s)
  
  end subroutine test_density
  


end module dcdft
