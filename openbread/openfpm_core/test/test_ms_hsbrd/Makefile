include ~/example.mk

CC=mpic++

LDIR = 

OBJ_ENZ = main.o

all:  model

debug: CFLAGS +=  -ggdb
debug: all

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

model: $(OBJ_ENZ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun -np 3 ./model;

.PHONY: clean all run

clean:
	rm -f *.o *~ core model

