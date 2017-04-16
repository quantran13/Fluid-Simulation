SOURCES_FILES=fluid.cu graphic.cu utility.cu fluid_kernels.cu
MAIN_FILE=main.cpp
BUILD_DIR=build/
INCLUDE_DIR=includes/
SRC_DIR=src/

CC=nvcc
CFLAGS=-g -I$(INCLUDE_DIR)
LIBFLAGS=-lm -lGL -lglut

OBJECTS=$(addprefix $(BUILD_DIR),$(notdir $(MAIN_FILE:.cpp=.o)))
OBJECTS := $(OBJECTS) $(addprefix $(BUILD_DIR),$(notdir $(SOURCES_FILES:.cu=.o)))

all: main

main: $(OBJECTS)
	$(CC) $(OBJECTS) -o main $(LIBFLAGS)

$(BUILD_DIR)%.o: $(SRC_DIR)%.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

$(BUILD_DIR)%.o: $(SRC_DIR)%.cu
	$(CC) $(CFLAGS) $^ -c -o $@

clean:
	rm -Rf main build/*
