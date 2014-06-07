name	  := Efficiency
FC        := gfortran
FFLAGS	  := -O
BINFOLDER := bin
WORKSPACE := tmp


libsuf := .lib
libfdr := xtt-lib
libsrc  := elliptic_tools.f90\
	   field_tools.f90\
	   constants.f90
libobj := $(addprefix $(WORKSPACE)/,$(libsrc:.f90=$(libsuf)))
libsrc := $(addprefix $(libfdr)/,$(libsrc))

libfile := $(WORKSPACE)/lib$(name).a
libfile_dep := $(libfile)($(notdir $(libobj)))


suf := .o
fdr := src
src	:= diagnose-real_data.f90\
	   diagnose-dry_test.f90\
	   diagnose-whs.f90
obj	:= $(addprefix $(WORKSPACE)/,$(src:.f90=$(suf)))
src := $(addprefix $(fdr)/,$(src))

.PHONY : all
all: | mkdir $(libfile_dep) main

.PHONY : mkdir
mkdir:
	mkdir $(WORKSPACE)
	mkdir $(BINFOLDER)

$(libfile_dep): $(libobj)
	ar -crs $@ $^


.PHONY : main
main: $(obj)
	for filename in $^; do \
		echo $${filename}; \
		mv $${filename} $(BINFOLDER)/; \
	done




$(WORKSPACE)/%$(libsuf): $(libfdr)/%.f90
	$(FC) $(FFLAGES) -c -o $@ $<

$(WORKSPACE)/%$(suf): $(fdr)/%.f90
	$(FC) $(FFLAGES) -o $@ $< $(libfile)



.PHONY : clean

clean:
	-rm *.mod
	-rm -r $(BINFOLDER)
	-rm -r $(WORKSPACE)

.PHONY : help


help:
	@echo "Help of $(name):"
	@echo "make"
	@echo "make clean"
