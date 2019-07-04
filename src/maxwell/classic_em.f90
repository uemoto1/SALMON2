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
subroutine classic_em
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use structures,           only: s_fdtd_system, s_fdtd_field
  use salmon_maxwell,       only: ls_fdtd_work,init_maxwell,calc_maxwell,finalize_maxwell
  use misc_routines,        only: get_wtime
  implicit none
  type(s_fdtd_system) :: fs
  type(s_fdtd_field)  :: ff
  type(ls_fdtd_work)  :: fw
  real(8)             :: elapsed_time
  
  elapsed_time=get_wtime()
  call init_maxwell(fs,ff,fw)
  call calc_maxwell(fs,ff,fw)
  call finalize_maxwell(fs,ff,fw)
  elapsed_time=get_wtime()-elapsed_time
  if(comm_is_root(nproc_id_global)) then
    write(*,'(A,f16.8)') " elapsed time [s] = ", elapsed_time
  end if
  
end subroutine classic_em
