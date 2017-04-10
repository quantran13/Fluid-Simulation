CC=gcc
CFLAGS = -g -Wall

all: main

main: fluid.c main.c fluid.h
	$(CC) $(CFLAGS) fluid.c main.c -o main -lm

clean:
	rm -Rf main
