CC = gcc
CFLAGS = -DNDEBUG -flto -I myopt
WARN = -Wwrite-strings -Wconversion -Wshadow -Wparentheses -Wlogical-op -Wunused -Wmissing-prototypes -Wmissing-declarations -Wdeclaration-after-statement -W -Wall -Wextra
OPT = -O3
LIBFLAGS = -lm

EXE = rmsdscan

SRC = main.c myopt/myopt.c pdb.c fit.c strtools.c proginfo.c 
DEPS = bool.h  myopt/myopt.h pdb.h fit.h strtools.h version.h proginfo.h
OBJ = main.o myopt/myopt.o pdb.o fit.o strtools.o proginfo.o

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








