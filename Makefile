# Makefile for spherIC

# Executable

BASE    = spherIC
EXT     = 
EXE     = spheric

# Compiler stuff

CC	= gcc
CFLAGS	= -O2 -Wall
LIBS	= -lm

# Object definition

OBJ	= $(EXE).o functions.o routines.o io.o

# Rules

$(EXE):	$(OBJ) Makefile
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LIBS)

clean:
	-rm -f *.o *~ $(EXE)

tar:
	cd ..; tar cvf - $(BASE)$(EXT)/*.c $(BASE)$(EXT)/*.h $(BASE)$(EXT)/Makefile $(BASE)$(EXT)/doc  > $(BASE)$(EXT).tar

# Dependencies

spheric.o: definitions.h functions.h routines.h io.h
functions.o: definitions.h functions.h routines.h
routines.o: definitions.h functions.h routines.h
io.o: io.h
