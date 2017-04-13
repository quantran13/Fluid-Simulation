CC=gcc
CFLAGS=-g -Wall -std=c99 -Iincludes/

all: main

main: src/fluid.c includes/fluid.h src/main.c src/graphic.c includes/graphic.h
	$(CC) $(CFLAGS) src/fluid.c src/main.c src/graphic.c -o main -lm -lGL -lglut

clean:
	rm -Rf main
