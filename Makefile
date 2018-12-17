utils=const.f90 random.f90 linalgebra.f90
src=${utils} PES.f90 DVR.f90
default:
	gfortran ${src} -llapack -lblas -o DVR
clean:
	rm *.mod
