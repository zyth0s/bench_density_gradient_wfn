
FFLAGS=-O3 -march=native -funroll-loops

libeval_gradrho.so: eval_gradrho.f90
		gfortran -Wall -shared -fPIC ${FFLAGS} -o libeval_gradrho.so eval_gradrho.f90
