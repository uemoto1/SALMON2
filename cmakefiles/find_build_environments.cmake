set(TARGET_SUFFIX ".cpu")

check_mpi_compiler(${CMAKE_Fortran_COMPILER} IS_MPI_COMPILER)
if (${IS_MPI_COMPILER})
  set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
endif ()

check_mpi_compiler(${CMAKE_C_COMPILER} IS_MPI_COMPILER)
if (${IS_MPI_COMPILER})
  set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
endif ()

if (USE_MPI)
  find_package(MPI REQUIRED)

  if (NOT DEFINED MPI_Fortran_FOUND)
    message(FATAL_ERROR "MPI Fortran compilers not found.")
  endif()

  if (NOT DEFINED MPI_C_FOUND)
    message(FATAL_ERROR "MPI C compilers not found.")
  endif()

  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

endif ()


set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(Fortran_FLAGS_General       "-cpp ${MPI_Fortran_COMPILE_FLAGS}")

set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-O3")
set(C_FLAGS_General             "${MPI_C_COMPILE_FLAGS}")

find_package(OpenMP REQUIRED)
set(OPENMP_FLAGS ${OpenMP_C_FLAGS})
