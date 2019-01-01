HERE= `pwd`
EXE2d=mc_VO2_2d
EXE3d=mc_VO2_3d

OBJS= 


##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)

FFLAG = -O3 -faggressive-loop-optimizations -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -Wall -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none


all: all compile
debug: debug compile

all: 

debug: FFLAG=$(DFLAG)


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(INCARGS) $(OBJS) $(EXE2d).f90 -o ~/.bin/$(EXE2d) $(LIBARGS)
	$(FC) $(FFLAG) $(INCARGS) $(OBJS) $(EXE3d).f90 -o ~/.bin/$(EXE3d) $(LIBARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" ~/.bin/$(EXE2d) "and" ~/.bin/$(EXE3d)

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@echo ""
