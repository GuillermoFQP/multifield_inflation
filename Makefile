# Configuration
FC      := gfortran
FLAGS   := -fopenmp -O3 -fdefault-real-8 -llapack -lblas
MODULES := multifield_globals.f90 multifield_utils.f90

# Programs
PROG1 := ps_full
PROG2 := ps_mode
PROG3 := efold_cmap
PROG4 := bkgd_trajs
PROG5 := potential_cmap
PROG6 := omega_eff
PROG7 := bkgd_onetraj
PROG8 := mode_injection
PROG9 := cov_mtx

# Default target
all: $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7) $(PROG8) $(PROG9)

# Build rules
$(PROG1): $(MODULES) $(PROG1).f90
	$(FC) $(FLAGS) $^ -o $@

$(PROG2): $(MODULES) $(PROG2).f90
	$(FC) $(FLAGS) $^ -o $@

$(PROG3): $(MODULES) $(PROG3).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG4): $(MODULES) $(PROG4).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG5): $(MODULES) $(PROG5).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG6): $(MODULES) $(PROG6).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG7): $(MODULES) $(PROG7).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG8): $(MODULES) $(PROG8).f90
	$(FC) $(FLAGS) $^ -o $@
	
$(PROG9): $(MODULES) $(PROG9).f90
	$(FC) $(FLAGS) $^ -o $@

# Clean targets
clean:
	@echo "Removing .o and .mod files..."
	@rm -f *.o *.mod

cleanout:
	@echo "Removing executables..."
	@rm -f $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7) $(PROG8) $(PROG9)

cleanall: clean cleanout

# Phony targets
.PHONY: all clean cleanout cleanall

