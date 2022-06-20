set -x

ifort -free -fpp -m64 -nothreads -assume byterecl -I/opt/intel/impi/2019.1.144/intel64/include -mkl=cluster \
-I$ADFHOME/build/libftl.build/intelmpi \
-I$ADFHOME/build/libscm_core.build/intelmpi -c read_adf.f90 

ifort -free -FR -fpp -m64 -nothreads -assume byterecl  -I/opt/intel/impi/2019.1.144/intel64/include -mkl=cluster \
-shared-intel -Bdynamic -lresolv  -lrfftw -lfftw \
-I$ADFHOME/src/lib/core/license -I$ADFHOME/src/include \
-L/opt/intel/mkl/lib/intel64 -Wl,--start-group -Bstatic -lmkl_intel_lp64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 \
-lmkl_sequential -lmkl_core -Wl,--end-group -lxcfun \
-Bdynamic -L$ADFBIN/TclTk/lib -ltcl8.6 \
-L/opt/intel/impi/2019.1.144/intel64/lib/release -lmpi \
-L/opt/intel/impi/2019.1.144/intel64/lib -lmpifort -lrt -lpthread -ldl  \
-o read_adf.x  read_adf.o            \
-L$ADFHOME/build/intelmpi/lib -lscm_core \
-L$ADFHOME/build/intelmpi/lib -lftl $ADFBIN/zlib/lib/libz.a -lresolv \
 $ADFHOME/bin.intelmpi/hwloc/lib/libhwloc.a

