CC=gcc
CFLAGS=-g -Wall -std=c99 -Iincludes/

all: main

main: src/fluid.c includes/fluid.h src/main.c src/graphic.c src/utility.c includes/utility.h includes/graphic.h
	$(CC) $(CFLAGS) src/fluid.c src/main.c src/graphic.c src/utility.c -o main -lm -lGL -lglut

clean:
	rm -Rf main
