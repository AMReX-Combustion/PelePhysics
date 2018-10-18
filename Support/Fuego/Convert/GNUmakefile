f77_sources := ckinterp.f cklib_d.f cktime.f xerror.f tranfit_d.f
f90_sources = strUtils.f90

f77_objects := $(f77_sources:%.f=%.o)
f90_objects := $(f90_sources:%.f90=%.o)

vpath %.f .

all: convert.exe

convert.exe: ${f77_objects} ${f90_objects} convert.o
	gfortran ${f77_objects} ${f90_objects} convert.o -o convert.exe

tranfit_d.o: tranfit_d.f string_utility.mod

string_utility.mod: ${f90_objects}

%.o: %.f
	gfortran -ffixed-line-length-100 -c $<

%.o: %.f90
	gfortran -c $<

clean:
	\rm -rf ${f77_objects} ${f90_objects} string_utility.mod *~ convert.o convert.exe
