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
SUBROUTINE init_ps(alat,brl,matrix_A)
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
use scf_data
use allocate_psl_sub
implicit none
real(8),intent(in) :: alat(3,3),brl(3,3),matrix_A(3,3)

if(iSCFRT==1)then
  if(comm_is_root(nproc_id_global))then
    print *,"----------------------------------- init_ps"
  end if
end if

call storevpp

Mps_all=0
select case(iperiodic)
case(0)
  call calcJxyz_all
  call calcuV
  call calcVpsl
  allocate(ppg%ekr_uV(ppg%nps,ppg%nlma,1))
  ppg%ekr_uV(:,:,1) = cmplx(ppg%uV)
case(3)
  select case(iflag_hartree)
  case(2)
    call calcVpsl_periodic(matrix_A,brl)
  case(4)
    call calcVpsl_periodic_FFTE
  end select
  call calcJxyz_all_periodic(alat,matrix_A)
  call calcuV
end select

return

END SUBROUTINE init_ps

SUBROUTINE dealloc_init_ps(ppg,ppg_all,ppn)
  use structures, only: s_pp_grid, s_pp_nlcc
  use salmon_global
  implicit none
  type(s_pp_grid) :: ppg,ppg_all
  type(s_pp_nlcc) :: ppn

  deallocate(ppg%jxyz, ppg%jxx, ppg%jyy, ppg%jzz, ppg%rxyz)
  deallocate(ppg%lma_tbl, ppg%ia_tbl)
  deallocate(ppg%rinv_uvu,ppg%uv,ppg%duv)

  deallocate(ppg_all%jxyz,ppg_all%jxx,ppg_all%jyy,ppg_all%jzz,ppg_all%rxyz)
  deallocate(ppg_all%lma_tbl, ppg_all%ia_tbl)
  deallocate(ppg_all%rinv_uvu,ppg_all%uv,ppg_all%duv)

  deallocate(ppg%Vpsl_atom)
  if(allocated(ppg%ekr_uV)) deallocate(ppg%ekr_uV)

  deallocate(ppn%rho_nlcc, ppn%tau_nlcc)
END SUBROUTINE dealloc_init_ps
