LN="-L/usr/lib64/ -llapack -lblas -I/usr/include/ -lfftw3" #
#gfortran -O0 -shared -fPIC -o Fortlib.so Fortlib.f90
#gfortran -O3 -shared -fPIC -o Fortlib.so Fortlib.f90
#gfortran -mtune=native -march=native -O3 -shared -fPIC -o Fortlib.so Fortlib.f90
#FC="gfortran -O0 -shared -fPIC"
FC="gfortran -mtune=native -march=native -fopenmp -Ofast -shared -fPIC"
OBJ="Fortlib.so"
SRC="Fortlib.f90"
$FC -o $OBJ $SRC $LN


