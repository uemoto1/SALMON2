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

  subroutine init_dcdft(dc)
    use structures
    use salmon_global, only: num_fragment
    use parallelization
    use mpi ! --> wrapper
    implicit none
    type(s_dcdft),intent(inout) :: dc
    !
    integer :: nproc,myrank
    integer :: comm_F,nproc_F,myrank_F
    integer :: npg,i,j,k,m,ierr
    
    dc%icomm_tot = nproc_group_global
    myrank = nproc_id_global
    nproc = nproc_size_global
    
    dc%n_frag = num_fragment(1)*num_fragment(2)*num_fragment(3)
    
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
    
!!!!!!!!! test
write(*,*) "dc test 1:",myrank_F,nproc_F,dc%i_frag,myrank,nproc
    
    ! Override global variables
    nproc_group_global = comm_F
    nproc_id_global = myrank_F
    nproc_size_global = nproc_F
    
  end subroutine


end module dcdft
