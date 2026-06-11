PROG = broyden
PATH_TO_LIKWID = /usr/local

CC = gcc -I$(PATH_TO_LIKWID)/include -L$(PATH_TO_LIKWID)/lib
OBJS = main.o helper.o utils.o

CFLAGS = -O3 -march=native -mavx -fopt-info-vec -Wall -Wextra -Wpedantic -Werror -flto -fopenmp
LFLAGS =-lm -llikwid

DISTFILES = *.c *.h Makefile resultados
DISTDIR = `basename ${PWD}`

.PHONY: all clean limpa benchmark

main.o: main.c ../helper.h ../utils.h
	$(CC) $(CFLAGS) -c $< -o $@
helper.o: ../helper.c ../helper.h
	$(CC) $(CFLAGS) -c $< -o $@
utils.o: ../utils.c ../utils.h
	$(CC) $(CFLAGS) -c $< -o $@
	
all: clean $(PROG)

benchmark: CC += -DBENCHMARKING -DLIKWID_PERFMON
benchmark: clean $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean limpa:
	@echo "Limpando ...."
	@rm -f *~ *.bak *.tmp
	@rm -f  $(PROG) *.o core a.out
