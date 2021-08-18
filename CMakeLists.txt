cmake_minimum_required(VERSION 3.20)

project("fatring" LANGUAGES "Fortran")

add_executable("fatring" "composition.F90" "constants.F90" "fatring.F90" "primes.F90")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -funroll-loops -fno-stack-arrays")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-inline-small-functions -fno-tree-pre")  # known to produce wrong code with gfortran 10.x
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -std=gnu -fimplicit-none -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -W -Wsurprising -Wcharacter-truncation -Wno-unused-function -fmodule-private -pedantic")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -ggdb -fcheck=all")
