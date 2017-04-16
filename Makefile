CC=nvcc
CFLAGS=-g -Iincludes/

all: main

main: src/fluid.cu includes/fluid.h src/main.cu src/graphic.cu src/utility.cu includes/utility.h includes/graphic.h
	$(CC) $(CFLAGS) src/fluid.cu src/main.cu src/graphic.cu src/utility.cu -o main -lm -lGL -lglut

clean:
	rm -Rf main
