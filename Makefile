CC=gcc
CFLAGS = -g -Wall -std=c99

all: main

main: fluid.c main.c fluid.h
	$(CC) $(CFLAGS) fluid.c main.c -o main -lm

clean:
	rm -Rf main
