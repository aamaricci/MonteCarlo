HERE= `pwd`
EXE2d=mc_VO2_2d
EXE3d=mc_VO2_3d

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
	@echo "Compile: 2D mc_SINGLE_SITE:"
	$(FC) $(FPPFLAG) -D_SITE $(FFLAG) $(INCARGS) $(OBJS) $(EXE2d).f90 -o ~/.bin/$(EXE2d) $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE2d)
	@echo ""
	@echo "Compile: 3D mc_SINGLE_SITE:"
	$(FC) $(FPPFLAG) -D_SITE $(FFLAG) $(INCARGS) $(OBJS) $(EXE3d).f90 -o ~/.bin/$(EXE3d) $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE3d)
	@echo ""
	@echo " ...................... done .............................. "

lattice: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: 2D mc_LATTICE:"
	$(FC) $(FPPFLAG) -D_ $(FFLAG) $(INCARGS) $(OBJS) $(EXE2d).f90 -o ~/.bin/$(EXE2d)_mclat $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE2d)_mclat
	@echo ""
	@echo "Compile: 3D mc_LATTICE:"
	$(FC) $(FPPFLAG) -D_ $(FFLAG) $(INCARGS) $(OBJS) $(EXE3d).f90 -o ~/.bin/$(EXE3d)_mclat $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE3d)_mclat
	@echo ""
	@echo " ...................... done .............................. "


debug: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: 2D mc_SINGLE_SITE DEBUG:"
	$(FC) $(FPPFLAG) -D_SITE $(DFLAG) $(INCARGS) $(OBJS) $(EXE2d).f90 -o ~/.bin/$(EXE2d)_debug $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE2d)_debug
	@echo ""
	@echo "Compile: 3D mc_SINGLE_SITE DEBUG:"
	$(FC) $(FPPFLAG) -D_SITE $(DFLAG) $(INCARGS) $(OBJS) $(EXE3d).f90 -o ~/.bin/$(EXE3d)_debug $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE3d)_debug
	@echo ""
	@echo " ...................... done .............................. "


lattice_debug: $(OBJS)
	@echo " ..................... compile ........................... "
	@echo ""
	@echo "Compile: 2D mc_LATTICE DEBUG:"
	$(FC) $(FPPFLAG) -D_ $(DFLAG) $(INCARGS) $(OBJS) $(EXE2d).f90 -o ~/.bin/$(EXE2d)_mclat_debug $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE2d)_mclat_debug
	@echo ""
	@echo "Compile: 3D mc_LATTICE DEBUG:"
	$(FC) $(FPPFLAG) -D_ $(DFLAG) $(INCARGS) $(OBJS) $(EXE3d).f90 -o ~/.bin/$(EXE3d)_mclat_debug $(LIBARGS)
	@echo ""
	@echo "created" ~/.bin/$(EXE3d)_mclat_debug
	@echo ""
	@echo " ...................... done .............................. "

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@echo ""
