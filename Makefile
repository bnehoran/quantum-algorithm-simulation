# adapted from http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
CC=gcc
WARNINGS=-Wall -Wextra
WARN_IGN=-Wno-unused-value -Wno-unused-variable \
-Wno-unused-parameter -Wno-unused-but-set-variable
CFLAGS=-std=c99 -g
ODIR=obj
LIBS=-lm -lgmp

DEPS = vector.h state.h operator.h pair.h register.h inverse.h bignum.h

_OBJ = vector.o state.o operator.o inverse.o register.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

default: run_tests
	run_tests

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(W)

run_tests: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core
	rm -f .err .error
	rm -f run_tests scratch

w warn: W = $(WARNINGS) $(WARN_IGN)
warnings: W = $(WARNINGS)
w warn warnings: clean default

vector_tests: vector.c vector.h pair.h bignum.h
	$(CC) -o vector_tests vector.c $(CFLAGS) $(LIBS) -nostartfiles -Wl,-eVector_main

s: scratch
	scratch

scratch: scratch.c
	gcc -o scratch scratch.c $(CFLAGS) $(LIBS)

# debug: run_testsD

# run_testsD: $(OBJ)
# 	gcc -o $@ $^ $(CFLAGS) $(LIBS) -g