CC = gcc
CFLAGS = -DNDEBUG -flto 
WARN = -Wwrite-strings -Wconversion -Wshadow -Wparentheses -Wlogical-op -Wunused -Wmissing-prototypes -Wmissing-declarations -Wdeclaration-after-statement -W -Wall -Wextra
OPT = -O3
LIBFLAGS = -lm

EXE = rmsdscan

SRC = main.c pdb.c fit.c strtools.c 
DEPS = bool.h  pdb.h fit.h strtools.h version.h
OBJ = main.o pdb.o fit.o strtools.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

rmsdscan: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(WARN) $(OPT) $(LIBFLAGS)

all:
	$(CC) $(CFLAGS) $(WARN) $(OPT) -o $(EXE) $(SRC) $(LIBFLAGS)

install:
	cp $(EXE) /usr/local/bin/$(EXE)

clean:
	rm -f *.o *~ myopt/*.o 








