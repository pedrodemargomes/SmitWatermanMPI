CFLAGS=-Wall -O3
CC=gcc

all: swparalelo swsequencial

clean:
	rm -f *.o swparalelo swsequencial

swparalelo: swParalelo.c
	$(CC) $< $(CFLAGS) -fopenmp -lm -o $@

swsequencial: swSequencial.c
	$(CC) $< $(CFLAGS) -o $@

debug: swParalelo.c swSequencial.c 
	$(CC) swParalelo.c  $(CFLAGS) -DDEBUG -fopenmp -lm -o swparalelodebug
	$(CC) swSequencial.c $(CFLAGS) -DDEBUG -o swsequencialdebug
