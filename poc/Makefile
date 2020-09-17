CC = gcc
FLAGS = -g -fopenmp
DEFS = -DTHREADS=8
LIBS = -lm -lrt -lopencore-amrnb

all: build

build: build.c
	${CC} ${FLAGS} ${DEFS} -o build build.c ${LIBS}

clean:
	rm -rf build
