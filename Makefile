HERE   = `pwd`
EXE=testISING_2D
DIREXE   = $(HERE)/test

OBJS     = 
OBJS_DEB =


##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)

FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none


all: all compile
debug: debug compile

all: 

debug: FFLAG=$(DFLAG)


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(INCARGS) $(OBJS) $(EXE).f90 -o test/$(EXE) $(LIBARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" test/$(EXE)

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@rm -f ${DIREXE}/${EXE}
	@echo ""
