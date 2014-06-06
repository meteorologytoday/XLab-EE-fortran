name	  := Efficiency
FC        := gfortran
FFLAGS	  := -O
BINFOLDER := bin
WORKSPACE := tmp

lib	:= lib$(name).a
libpath := xtt-lib
libsrc  := elliptic_tools.f90\
	   field_tools.f90\
	   constants.f90

libobj	:= $(addprefix $(WORKSPACE)/,$(libsrc:.f90=.lib.out))
libsrc	:= $(addprefix $(libpath)/,$(libsrc))

srcpath := src
src	:= diagnose-real_data.f90\
	   diagnose-dry_test.f90\
	   diagnose-whs.f90

obj	:= $(addprefix $(WORKSPACE)/,$(src:.f90=.out))
src := $(addprefix $(srcpath)/,$(src))



.PHONY : all

all: mkdir $(libobj) $(obj)
	mv *.out $(BINFOLDER)
	mv $(lib) $(BINFOLDER)

mkdir:
	mkdir $(WORKSPACE)
	mkdir $(BINFOLDER)

$(WORKSPACE)/%.lib.out : $(libpath)/%.f90
	@echo "Making making " $@
	$(FC) $(FFLAGS) -c $@

$(WORKSPACE)/%.out : $(srcpath)/%.f90 $(libobj)
	@echo "Now making " $@
	$(FC) $(FFLAGS) -c -o $@ $?


.PHONY : clean

clean:
	-rm -f *.o
	-rm -f *.out
	-rm -f *.mod
	-rm -f $(lib)
	-rm -r $(BINFOLDER)
	-rm -r $(WORKSPACE)

.PHONY : help


help:
	@echo "Help of $(name):"
	@echo "make"
	@echo "make clean"
