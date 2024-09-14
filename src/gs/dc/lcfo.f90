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
module lcfo ! DC-LCFO: Phys. Rev. B 95, 045106 (2017).
  implicit none
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
      integer :: id_src,id_dst,ifrag_src,ifrag_dst,dvec(3),length(3),dsp_send(3),dsp_recv(3)
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
    n_mat = dc%nstate_frag * dc%n_frag
    
    call init_lcfo
    
    call calc_basis
    
    call hpsi_basis
    
    call calc_hamiltonian_matrix
    
    do ispin=1,nspin
      call eigen_dsyev(mat_H(1:n_mat(ispin),1:n_mat(ispin),ispin),esp_tot(1:n_mat(ispin),ispin) &
      &               ,mat_V(1:n_mat(ispin),1:n_mat(ispin),ispin))
    end do
  

if(dc%id_tot==0) then!!!!!!!!! test_dcdft
  do i=1,n_mat(1)
    write(777,*) i,esp_tot(i,1)
  end do
end if!!!!!!!!! test_dcdft
    
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
            halo(i)%id_dst = id_array(ifrag)
            halo(i)%ifrag_dst = ifrag
          end if
        ! src neighbor (-)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) - halo(i)%dvec(1:3)*dc%nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_src < 0) then
            halo(i)%id_src = id_array(ifrag)
            halo(i)%ifrag_src = ifrag
          end if
        end do ! ifrag
        if(halo(i)%id_dst < 0 .or. halo(i)%id_src < 0) stop "DC-LCFO: dst, src"
      end do
      end do
      end do
      n_halo = i
      
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
      
    ! f_basis <-- | \bar{\phi} >
      wrk_array = 0d0
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if( ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. iz <= dc%nxyz_domain(3) &
        & .and. energy%esp(io,1,ispin) - system%mu < energy_cut ) then
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
          if(lambda(io,ispin) > lambda_cut) then
            i = i + 1 ! count # of basis functions
            do jo=1,dc%nstate_frag
              f_basis(:,:,:,ispin,i) = f_basis(:,:,:,ispin,i) &
              & + wrk_array(:,:,:,ispin,jo) * mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
            end do
          end if
        end do ! io
        nb(ispin) = i ! # of basis functions
      end do ! ispin
      
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
      
    ! hat_H <-- Hamiltonian matrix < lambda | H | lambda >
      nmax = maxval(n_mat)
      allocate(mat_H(nmax,nmax,nspin),mat_V(nmax,nmax,nspin),esp_tot(nmax,nspin))
      mat_V = 0d0
      if(dc%id_frag==0) then
        ifrag = dc%i_frag
        l = dc%nxyz_domain
        do ispin=1,nspin
        ! diagonal block
          do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
          do jo=1,n_basis(ifrag,ispin) ; j = index_basis(jo,ifrag,ispin)
            mat_V(i,j,ispin) = mat_V(i,j,ispin) &
            & + sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
          end do
          end do
        ! off-diagonal block
          do i_halo=1,n_halo ; jfrag = halo(i_halo)%ifrag_src ! src fragment (recv)
            do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
            do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
            ! mat_H_local(jo,io) == < lambda_{jfrag,jo} | H | lambda_{ifrag,io} >
              wrk = 0.5d0* halo(i_halo)%mat_H_local(jo,io,ispin)
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
  
  end subroutine dc_lcfo

end module lcfo
