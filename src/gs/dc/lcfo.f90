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

! DC-LCFO method [Phys. Rev. B 95, 045106 (2017).]

module lcfo
  implicit none
  
  private
  public :: dc_lcfo, restart_rt_from_data_dcdft
  
  character(32),parameter :: binfile_wf = "wavefunctions.bin", &
  &                          binfile_rg = "rgrid_index.bin", &
  &                          binfile_bf = "basis_functions.bin", &
  &                          binfile_hl = "hamiltonian_local.bin"
  
contains

  subroutine dc_lcfo(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_summation
    use eigen_subdiag_sub, only: eigen_dsyev
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),      intent(in) :: stencil
    type(s_pp_grid),      intent(in) :: ppg
    type(s_dft_energy),   intent(in) :: energy
    type(s_scalar),       intent(in) :: V_local(system%nspin)
    type(s_orbital),      intent(in) :: spsi
    type(s_orbital)                  :: shpsi,sttpsi
    type(s_sendrecv_grid)            :: srg
    type(s_dcdft)                    :: dc
    !
    type halo_info
      integer :: id_src,id_dst,ifrag_src,dvec(3),length(3),dsp_send(3),dsp_recv(3)
      real(8),allocatable :: buf_send(:,:,:,:,:),buf_recv(:,:,:,:,:),mat_H_local(:,:,:)
    end type halo_info
    !
    type(halo_info) :: halo(26) ! 26 = 3^3-1
    integer :: nspin,n_halo
    integer :: n_basis(dc%n_frag,system%nspin), n_mat(system%nspin)
    integer :: index_basis(dc%nstate_frag,dc%n_frag,system%nspin)
    real(8) :: hvol
    real(8),allocatable :: f_basis(:,:,:,:,:),hf(:,:,:,:,:),wrk_array(:,:,:,:,:), &
    & mat_H(:,:,:),mat_V(:,:,:),esp_tot(:,:)
    !
    integer :: i,j,n,ix,iy,iz,io,jo,ispin,ifrag,jfrag

integer :: nnn !!!!!!!!! test_lcfo
    
    hvol = system%hvol
    nspin = system%nspin
    
    call init_lcfo
    
    call calc_basis
    
    call hpsi_basis
    
    call calc_hamiltonian_matrix
    
    do ispin=1,nspin
      call eigen_dsyev(mat_H(1:n_mat(ispin),1:n_mat(ispin),ispin),esp_tot(1:n_mat(ispin),ispin) &
      &               ,mat_V(1:n_mat(ispin),1:n_mat(ispin),ispin))
    end do
  
    call output

    deallocate(f_basis)
    do i=1,n_halo
      if(allocated(halo(i)%mat_H_local)) deallocate(halo(i)%mat_H_local)
    end do
    deallocate(mat_H,mat_V,esp_tot)
    
  contains
  
    subroutine init_lcfo
      use salmon_global, only: num_fragment
      implicit none
      integer :: lx,ly,lz
      integer,dimension(3) :: nh,ir1,ir2,d
      integer,dimension(dc%n_frag) :: id_array, id_tmp
      
      id_tmp = 0
      if(dc%id_frag==0) id_tmp(dc%i_frag) = dc%id_tot + 1
      call comm_summation(id_tmp,id_array,dc%n_frag,dc%icomm_tot)
      id_array = id_array - 1
      
      nh = 0
      do n=1,3 ! x,y,z
        if(dc%nxyz_buffer(n) > dc%nxyz_domain(n)) stop "DC-LCFO: buffer > domain"
        if(num_fragment(n) > 1) nh(n) = 1
      end do
      
      i = 0
      do lx=-nh(1),nh(1)
      do ly=-nh(2),nh(2)
      do lz=-nh(3),nh(3)
        if(lx==0 .and. ly==0 .and. lz==0) cycle
        i = i + 1
        halo(i)%dvec(1:3) = [lx, ly, lz]
        halo(i)%id_dst = -1
        halo(i)%id_src = -1
        do ifrag=1,dc%n_frag
        ! dc%ixyz_frag: r-grid index of the fragment origin
          ir1(1:3) = dc%ixyz_frag(1:3,ifrag) ! position of fragment ifrag
        ! dst neighbor (+)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) + halo(i)%dvec(1:3)*dc%nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_dst < 0) then
            halo(i)%id_dst = id_array(ifrag) ! process ID of the communication destination
          end if
        ! src neighbor (-)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) - halo(i)%dvec(1:3)*dc%nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_src < 0) then
            halo(i)%id_src = id_array(ifrag) ! process ID of the communication source
            halo(i)%ifrag_src = ifrag
          end if
        end do ! ifrag
        if(halo(i)%id_dst < 0 .or. halo(i)%id_src < 0) stop "DC-LCFO: dst, src"
      end do
      end do
      end do
      n_halo = i ! # of the halo regions (neighbor fragments)
      
      do i=1,n_halo
        do n=1,3 ! x,y,z
          select case (halo(i)%dvec(n))
          case(0)
            halo(i)%length(n) = dc%nxyz_domain(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = 0
          case(1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = dc%nxyz_domain(n) - dc%nxyz_buffer(n)
            halo(i)%dsp_recv(n) = dc%nxyz_domain(n) + dc%nxyz_buffer(n)
          case(-1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = dc%nxyz_domain(n)
          end select
        end do
      end do
    
    end subroutine init_lcfo
    
    subroutine calc_basis
      use salmon_global, only: energy_cut,lambda_cut
      implicit none
      integer :: nb(nspin),itmp(dc%n_frag,nspin)
      real(8),dimension(dc%nstate_frag,dc%nstate_frag,system%nspin) :: mat_S,mat_U
      real(8),dimension(dc%nstate_frag,system%nspin) :: lambda
      
      allocate(f_basis  (dc%nxyz_domain(1),dc%nxyz_domain(2),dc%nxyz_domain(3),nspin,dc%nstate_frag))
      allocate(wrk_array(dc%nxyz_domain(1),dc%nxyz_domain(2),dc%nxyz_domain(3),nspin,dc%nstate_frag))
      
    ! f_basis <-- | \bar{\phi} > (projected fragment orbitals)
      wrk_array = 0d0
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if( ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. iz <= dc%nxyz_domain(3) &
        & .and. energy%esp(io,1,ispin) - system%mu < energy_cut ) then ! energy cutoff
          wrk_array(ix,iy,iz,ispin,io) = spsi%rwf(ix,iy,iz,ispin,io,1,1) ! | \phi > @ core domain
        end if
      end do
      end do
      end do
      end do
      end do
      call comm_summation(wrk_array,f_basis,product(dc%nxyz_domain)*nspin*dc%nstate_frag,info%icomm_rko)
      
    ! mat_S <-- S_{ij} = < \bar{\phi}_i | \bar{\phi}_j > (overlap matrix)
      do ispin=1,nspin
      do io=1,dc%nstate_frag
      do jo=1,dc%nstate_frag
        mat_S(io,jo,ispin) = sum(f_basis(:,:,:,ispin,io)*f_basis(:,:,:,ispin,jo)) * hvol
      end do
      end do
      end do
      
    ! diagonalize mat_S
      do ispin=1,nspin
        call eigen_dsyev(mat_S(:,:,ispin),lambda(:,ispin),mat_U(:,:,ispin))
      end do
      
    ! f_basis <-- | lambda > (basis functions)
      wrk_array = f_basis
      f_basis = 0d0
      do ispin=1,nspin
        i = 0 ! count # of basis functions
        do io=dc%nstate_frag,1,-1
          if( lambda(io,ispin) > lambda_cut ) then ! cutoff for the eigenvalues of the overlap matrix
            i = i + 1 ! count # of basis functions
            do jo=1,dc%nstate_frag
              f_basis(:,:,:,ispin,i) = f_basis(:,:,:,ispin,i) &
              & + wrk_array(:,:,:,ispin,jo) * mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
            end do
          end if
        end do ! io
        nb(ispin) = i ! # of basis functions
      end do ! ispin
      
    ! Gram–Schmidt orthonormalization
      wrk_array = f_basis
      do ispin=1,nspin
        do io=1,nb(ispin)
          do jo=1,io-1
            wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
            & - f_basis(:,:,:,ispin,jo) * sum(f_basis(:,:,:,ispin,jo)*wrk_array(:,:,:,ispin,io)) &
            & / sum(f_basis(:,:,:,ispin,jo)*f_basis(:,:,:,ispin,jo))
          end do
          wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
          & / sqrt( sum(wrk_array(:,:,:,ispin,io)*wrk_array(:,:,:,ispin,io)) * hvol )
        end do
      end do ! ispin
      f_basis = wrk_array
      
    ! sttpsi <-- f_basis == | lambda > (basis functions)
      sttpsi%rwf = 0d0
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if( ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. iz <= dc%nxyz_domain(3) ) then
          sttpsi%rwf(ix,iy,iz,ispin,io,1,1) = f_basis(ix,iy,iz,ispin,io)
        end if
      end do
      end do
      end do
      end do
      end do
      
    ! n_basis: # of basis functions
      itmp = 0
      if(dc%id_frag==0) itmp(dc%i_frag,1:nspin) = nb(1:nspin)
      call comm_summation(itmp,n_basis,dc%n_frag*nspin,dc%icomm_tot)
      index_basis = 0
      do ispin=1,nspin
        i = 0
        do ifrag=1,dc%n_frag
          do io=1,n_basis(ifrag,ispin)
            i = i + 1
            index_basis(io,ifrag,ispin) = i ! index_basis: index for the total matrix
          end do
        end do
        n_mat(ispin) = i ! n_mat: dimension of the total matrix
      end do

if(dc%id_frag==0) then !!!!!!!!! test_lcfo
write(*,*) "test_lcfo n_basis", nb
nnn=4000
nnn=nnn+dc%i_frag
write(nnn,*) "# lambda"
do io=dc%nstate_frag,1,-1
  write(nnn,*) lambda(io,1)
end do
nnn=1000
nnn=nnn+dc%i_frag
write(nnn,*) "# i, j, < lambda_i | lambda_j >"
do io=1,dc%nstate_frag
do jo=1,dc%nstate_frag
  write(nnn,*) io,jo, sum(f_basis(:,:,:,1,io)*f_basis(:,:,:,1,jo)) * hvol
end do
end do
end if !!!!!!!!! test_lcfo
      
      deallocate(wrk_array)
      
    end subroutine calc_basis
    
    subroutine hpsi_basis
      use hamiltonian, only: hpsi
      implicit none
      
      allocate(hf       (lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))
      allocate(wrk_array(lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))
      
    ! shpsi <-- H | lambda > (Hamiltonian operation)
      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)
      
    ! hf <-- shpsi == H | lambda >
      wrk_array = 0d0
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk_array(ix,iy,iz,ispin,io) = shpsi%rwf(ix,iy,iz,ispin,io,1,1)
      end do
      end do
      end do
      end do
      end do
      call comm_summation(wrk_array,hf,product(lg%num)*nspin*dc%nstate_frag,info%icomm_rko)
      
      deallocate(wrk_array)
      
    end subroutine hpsi_basis
    
    subroutine calc_hamiltonian_matrix
      use communication, only: comm_isend, comm_irecv, comm_wait_all
      implicit none
      integer :: i_halo,d(3),l(3),nmax
      integer :: itag_send,itag_recv
      integer,dimension(n_halo) :: ireq_send,ireq_recv
      real(8) :: wrk
            
      if(dc%id_frag==0) then
      ! halo communication for f_basis == | lambda >
        do i_halo=1,n_halo
          l = halo(i_halo)%length
          d = halo(i_halo)%dsp_send
          allocate(halo(i_halo)%buf_send(l(1),l(2),l(3),nspin,dc%nstate_frag))
          allocate(halo(i_halo)%buf_recv(l(1),l(2),l(3),nspin,dc%nstate_frag))
          do iz=1,l(3)
          do iy=1,l(2)
          do ix=1,l(1)
            halo(i_halo)%buf_send(ix,iy,iz,:,:) = f_basis(d(1)+ix,d(2)+iy,d(3)+iz,:,:)
          end do
          end do
          end do
write(*,*) "test_lcfo: sendrecv",i_halo,dc%i_frag,dc%id_tot,halo(i_halo)%id_dst,halo(i_halo)%id_src !!!!!!!!! test_lcfo
        ! MPI_ISEND
          itag_send = dc%i_frag
          ireq_send(i_halo) = comm_isend(halo(i_halo)%buf_send,halo(i_halo)%id_dst,itag_send,dc%icomm_tot)
        ! MPI_IRECV
          itag_recv = halo(i_halo)%ifrag_src
          ireq_recv(i_halo) = comm_irecv(halo(i_halo)%buf_recv,halo(i_halo)%id_src,itag_recv,dc%icomm_tot)
        end do ! i_halo=1,n_halo
        call comm_wait_all(ireq_recv)
        call comm_wait_all(ireq_send)
      ! halo(i_halo)%mat_H_local: local Hamiltonian matrix < lambda | H | lambda > (off-diagonal block)
        do i_halo=1,n_halo
          l = halo(i_halo)%length
          d = halo(i_halo)%dsp_recv
          allocate(halo(i_halo)%mat_H_local(dc%nstate_frag,dc%nstate_frag,system%nspin))
          halo(i_halo)%mat_H_local = 0d0
          do ispin=1,nspin
          do io=1,dc%nstate_frag
          do jo=1,dc%nstate_frag
            do iz=1,l(3)
            do iy=1,l(2)
            do ix=1,l(1)
              halo(i_halo)%mat_H_local(io,jo,ispin) = halo(i_halo)%mat_H_local(io,jo,ispin) &
              & + halo(i_halo)%buf_recv(ix,iy,iz,ispin,io) * hf(d(1)+ix,d(2)+iy,d(3)+iz,ispin,jo) * hvol
            end do
            end do
            end do
          end do
          end do
          end do
          deallocate(halo(i_halo)%buf_send,halo(i_halo)%buf_recv)
        end do
      end if ! dc%id_frag==0
      
if(dc%id_frag==0) then !!!!!!!!! test_lcfo
nnn=2000
nnn=nnn+dc%i_frag
l = dc%nxyz_domain
write(nnn,*) "# i,j,mat_H(i,j)"
do io=1,dc%nstate_frag
do jo=1,dc%nstate_frag
  write(nnn,*) io,jo, sum(f_basis(1:l(1),1:l(2),1:l(3),1,io)*hf(1:l(1),1:l(2),1:l(3),1,jo)) * hvol
end do
end do
nnn=3000
nnn=nnn+dc%i_frag
l = dc%nxyz_domain
write(nnn,*) "# i_halo,i,j,mat_H(i,j) (inter-fragment)"
do i_halo=1,n_halo
do io=1,dc%nstate_frag
do jo=1,dc%nstate_frag
  write(nnn,*) i_halo,io,jo,halo(i_halo)%mat_H_local(io,jo,1)
end do
end do
end do
end if !!!!!!!!! test_lcfo
      
    ! mat_H <-- total Hamiltonian matrix < lambda | H | lambda >
      nmax = maxval(n_mat)
      allocate(mat_H(nmax,nmax,nspin),mat_V(nmax,nmax,nspin),esp_tot(nmax,nspin))
      mat_V = 0d0
      if(dc%id_frag==0) then
        ifrag = dc%i_frag
        l = dc%nxyz_domain
        do ispin=1,nspin
        ! diagonal block < lambda_{ifrag,io} | H | lambda_{ifrag,jo} >
          do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
          do jo=1,n_basis(ifrag,ispin) ; j = index_basis(jo,ifrag,ispin)
            mat_V(i,j,ispin) = mat_V(i,j,ispin) &
            & + sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
          end do
          end do
        ! off-diagonal block < lambda_{jfrag,jo} | H | lambda_{ifrag,io} >
          do i_halo=1,n_halo ; jfrag = halo(i_halo)%ifrag_src ! src fragment (recv)
            do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
            do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
            ! mat_H_local(jo,io) == < lambda_{jfrag,jo} | H | lambda_{ifrag,io} >
              wrk = 0.5d0* halo(i_halo)%mat_H_local(jo,io,ispin) ! 0.5d0* : for double counting
              mat_V(j,i,ispin) = mat_V(j,i,ispin) + wrk
              mat_V(i,j,ispin) = mat_V(i,j,ispin) + wrk
            end do
            end do
          end do
        end do ! ispin=1,nspin
      end if ! dc%id_frag==0
      call comm_summation(mat_V,mat_H,nmax*nmax*nspin,dc%icomm_tot)
      deallocate(hf)
      
    end subroutine calc_hamiltonian_matrix
    
    subroutine output
      use salmon_global, only: base_directory, sysname, unit_energy
      use filesystem, only: get_filehandle
      use inputoutput, only: uenergy_from_au
      implicit none
      integer :: iunit,i_halo
      character(256) :: filename
      
    ! total system data
      if(dc%id_tot==0) then
      ! eigen.data
        iunit = get_filehandle()
        filename = trim(dc%base_directory)//trim(sysname)//"_eigen.data" ! @ ./data_dcdft/total/
        open(iunit,file=filename)
        write(iunit,'("#esp: single-particle energies (eigen energies) calculated by DC-LCFO method")')
        write(iunit,'("#io: orbital index")')
        select case(unit_energy)
        case('au','a.u.')
          write(iunit,'("# 1:io, 2:esp[a.u.]")')
        case('ev','eV')
          write(iunit,'("# 1:io, 2:esp[eV]")')
        end select
        do ispin=1,nspin
          write(iunit,'("# spin=",1x,i5)') ispin
          do i=1,n_mat(ispin)
            write(iunit,'(1x,i5,e26.16e3)') i,esp_tot(i,ispin)*uenergy_from_au
          end do
        end do
        close(iunit)
      ! coefficients of the wavefunctions
        iunit = get_filehandle()
        filename = trim(dc%base_directory)//binfile_wf ! @ ./data_dcdft/total/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) dc%n_frag, nspin, dc%nstate_frag
        write(iunit) n_mat(1:nspin)
        write(iunit) n_basis(1:dc%n_frag,1:nspin)
        write(iunit) index_basis(1:dc%nstate_frag,1:dc%n_frag,1:nspin)
        do ispin=1,nspin
          write(iunit) mat_V(1:n_mat(ispin),1:n_mat(ispin),ispin)
        end do
        close(iunit)
      end if
      
    ! fragment data
      if(dc%id_frag==0) then
      ! r-grid index
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_rg ! base_directory==./data_dcdft/fragments/dc%i_frag/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) lg%num(1:3), dc%lg_tot%num(1:3)
        do n=1,3 ! x,y,z
          write(iunit) dc%jxyz_tot(1:lg%num(n),n)
        end do
        close(iunit)
      ! basis functions | lambda >
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_bf ! base_directory==./data_dcdft/fragments/dc%i_frag/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) dc%nxyz_domain(1:3),nspin,dc%nstate_frag
        write(iunit) n_basis(dc%i_frag,1:nspin) ! # of basis functions
        write(iunit) f_basis(1:dc%nxyz_domain(1),1:dc%nxyz_domain(2),1:dc%nxyz_domain(3) &
        & ,1:nspin,1:dc%nstate_frag) ! basis functions | lambda >
        close(iunit)
      ! local hamiltonian matrix
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_hl
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
      end if
      
    end subroutine output
  
  end subroutine dc_lcfo

!===================================================================================================================================

! TDDFT & yn_dc==y : conventional TDDFT but wavefunctions are reconstructed from DC-LCFO data
  subroutine restart_rt_from_data_dcdft(lg,mg,system,info,spsi)
    use communication, only: comm_is_root, comm_summation, comm_bcast
    use filesystem, only: get_filehandle
    use salmon_global,only: num_fragment
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_orbital)                  :: spsi
    !
    character(32),parameter :: bdir_tot='./data_dcdft/total/',bdir_frag='./data_dcdft/fragments/'
    character(256) :: filename
    integer :: iunit, n_frag, nspin, nstate_frag, nmax
    integer :: i,j,jfrag,ispin,io,jo,ix,iy,iz,ix_tot,iy_tot,iz_tot,n
    integer,dimension(3) :: lgnum_frag,lgnum_tmp,nxyz_domain
    !
    integer,allocatable :: n_mat(:),n_basis(:,:),index_basis(:,:,:),jxyz_tot(:,:)
    real(8),allocatable :: f_basis(:,:,:,:,:),mat_V(:,:,:),wrk1(:,:,:),wrk2(:,:,:)
    
  ! read coefficients of the wavefunctions ./data_dcdft/total/wavefunctions.bin
    allocate(n_mat(system%nspin))
    if(comm_is_root(info%id_rko)) then
      write(*,*) "TDDFT & yn_dc==y : conventional TDDFT but wavefunctions are reconstructed from DC-LCFO data"
      write(*,*) "read from ./data_dcdft directory"
      iunit = get_filehandle()
      filename = trim(bdir_tot)//binfile_wf ! @ ./data_dcdft/total/
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) n_frag, nspin, nstate_frag
      if( n_frag /= product(num_fragment) .or. nspin /= system%nspin ) stop "data_dcdft: input mismatch"
      read(iunit) n_mat(1:nspin)
    end if
    nspin = system%nspin
    call comm_bcast(n_frag,info%icomm_rko)
    call comm_bcast(nstate_frag,info%icomm_rko)
    call comm_bcast(n_mat,info%icomm_rko)
    allocate(n_basis(n_frag,nspin))
    allocate(index_basis(nstate_frag,n_frag,nspin))
    nmax = maxval(n_mat)
    allocate(mat_V(nmax,nmax,nspin))
    if(comm_is_root(info%id_rko)) then
      read(iunit) n_basis(1:n_frag,1:nspin)
      read(iunit) index_basis(1:nstate_frag,1:n_frag,1:nspin)
      do ispin=1,nspin
        read(iunit) mat_V(1:n_mat(ispin),1:n_mat(ispin),ispin)
      end do
      close(iunit)
    end if
    call comm_bcast(n_basis,info%icomm_rko)
    do ispin=1,nspin
      call comm_bcast(index_basis(1:nstate_frag,1:n_frag,ispin),info%icomm_rko)
    end do
    call comm_bcast(mat_V,info%icomm_rko)
    
  ! read fragment basis functions ./data_dcdft/fragments/*/basis_functions.bin
    if(info%isize_rko < n_frag) stop "yn_dc=y: MPI size is too small."
    n = info%isize_rko / n_frag
    jfrag = -1
    do j=1,n_frag
      if( j*n == (info%id_rko+1) ) then
        jfrag = j
        exit
      end if
    end do
    if(jfrag > 0) then ! myrank (info%id_rko) <--> jfrag
    ! r-grid index
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_rg
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) lgnum_frag(1:3), lgnum_tmp(1:3)
      if( any( lgnum_tmp /= lg%num ) ) stop "data_dcdft: input mismatch (lg)"
      allocate(jxyz_tot(maxval(lgnum_frag),3))
      do n=1,3 ! x,y,z
        read(iunit) jxyz_tot(1:lgnum_frag(n),n)
      end do
      close(iunit)
    ! basis functions | lambda >
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_bf
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) nxyz_domain(1:3),i,j ! i,j: dummy
      read(iunit) lgnum_tmp(1:nspin) ! dummy
      if(i /= nspin .or. j /= nstate_frag .or. any( lgnum_tmp(1:nspin) /= n_basis(jfrag,1:nspin) ) ) then
        stop "data_dcdft: input mismatch (basis_functions.bin)"
      end if
      allocate(f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:nstate_frag))
      read(iunit) f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:nstate_frag)
      close(iunit)
    end if
    
  ! r-grid wavefunctions
    allocate(wrk1(lg%num(1),lg%num(2),lg%num(3)))
    allocate(wrk2(lg%num(1),lg%num(2),lg%num(3)))
    do ispin=1,nspin
    do io=1,system%no
      wrk1 = 0d0
      if(jfrag > 0) then ! myrank (info%id_rko) <--> jfrag
        do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
        do iz=1,nxyz_domain(3); iz_tot = jxyz_tot(iz,3)
        do iy=1,nxyz_domain(2); iy_tot = jxyz_tot(iy,2)
        do ix=1,nxyz_domain(1); ix_tot = jxyz_tot(ix,1)
          wrk1(ix_tot,iy_tot,iz_tot) = wrk1(ix_tot,iy_tot,iz_tot) &
          & + f_basis(ix,iy,iz,ispin,jo) * mat_V(j,io,ispin)
        end do
        end do
        end do
        end do
      end if
      call comm_summation(wrk1,wrk2,product(lg%num(1:3)),info%icomm_rko)
      if(info%io_s <= io .and. io <= info%io_e) then
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          spsi%zwf(ix,iy,iz,ispin,io,1,1) = wrk2(ix,iy,iz)
        end do
        end do
        end do        
      end if
    end do
    end do
    
    deallocate(n_mat,n_basis,index_basis,mat_V,wrk1,wrk2)
  end subroutine restart_rt_from_data_dcdft
  

end module lcfo
