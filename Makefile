HERE= `pwd`
EXE=mc_VO2_2d

OBJS= 

##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)

FFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -Wall -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
FPPFLAG =-cpp

all: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: mc_LATTICE:"
	$(FC) $(FPPFLAG) -D_ $(FFLAG) $(INCARGS) $(OBJS) $(EXE).f90 -o ~/.bin/$(EXE) $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE)
	@echo " ...................... done .............................. "

debug: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: mc_LATTICE DEBUG:"
	$(FC) $(FPPFLAG) -D_ $(DFLAG) $(INCARGS) $(OBJS) $(EXE).f90 -o ~/.bin/$(EXE)_debug $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE)_debug
	@echo " ...................... done .............................. "


site: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: mc_SINGLE_SITE:"
	$(FC) $(FPPFLAG) -D_SITE $(FFLAG) $(INCARGS) $(OBJS) $(EXE).f90 -o ~/.bin/$(EXE)_site $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE)_site
	@echo " ...................... done .............................. "




site_debug: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: mc_SINGLE_SITE DEBUG:"
	$(FC) $(FPPFLAG) -D_SITE $(DFLAG) $(INCARGS) $(OBJS) $(EXE).f90 -o ~/.bin/$(EXE)_site_debug $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE)_site_debug
	@echo " ...................... done .............................. "

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@echo ""
