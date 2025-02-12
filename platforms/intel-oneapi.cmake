### Intel Compiler for AVX-512 processors (since Skylake-SP)
set(ARCH_FLAGS                  "-xCORE-AVX512")
set(SIMD_SET                    "avx512")
set(OPENMP_FLAGS                "-qopenmp")
set(LAPACK_VENDOR_FLAGS         "-qmkl=parallel")
set(ScaLAPACK_VENDOR_FLAGS      "-qmkl=cluster")
set(Fortran_PP_FLAGS            "-fpp")

set(CMAKE_Fortran_COMPILER      "mpiifx")
set(CMAKE_C_COMPILER            "mpiicx")

set(General_Fortran_FLAGS       "-nogen-interface -std03 -warn all -diag-disable 5194,5462,6477,6717,6916,7025,7346,7373,7374,7416,7712,8430 -ansi-alias -fno-alias -qoverride-limits")
set(General_C_FLAGS             "-Wall -diag-disable=10388 -restrict -ansi-alias -fno-alias")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g ${General_C_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_C_FLAGS}")

set(USE_MPI_DEFAULT                        ON)
#set(USE_OPT_EXPLICIT_VECTORIZATION_DEFAULT ON)


########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Intel AVX-512 processors")
set(CMAKE_SYSTEM_PROCESSOR "skylake")
