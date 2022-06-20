ifort -fpp -qopenmp -DOMP -o write_wavet.x write_wavet.f90
