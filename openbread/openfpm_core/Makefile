include ~/example.mk

CC=mpic++
LDIR =
#OBJ_ENZ = pyhsbrd.o
PYTHON_INCLUDE = /ssoft/spack/cornalin/v1/opt/spack/linux-rhel7-x86_E5v4_Mellanox/gcc-5.4.0/python-3.6.1-moc7d7wihdugzp5fl4ew7za5muzwvrn6/include/python3.6m
#PYTHON_INCLUDE = /usr/local/include/python3.5m

CFLAGS += -fPIC

all: pyhsbrd

debug: CFLAGS += -ggdb
debug: all

pyhsbrd:
			swig -python -c++ pyhsbrd.i
			$(CC)  -O3 -g -c --std=c++11  pyhsbrd.cpp \
			$< $(INCLUDE_PATH) $(CFLAGS) $(LIBS_PATH) $(LIBS) \
			-DHAVE_CONFIG_H -I -I$(PYTHON_INCLUDE)

			$(CC)  -O3 -g -c --std=c++11 pyhsbrd_wrap.cxx \
			$< $(INCLUDE_PATH) $(CFLAGS) $(LIBS_PATH) $(LIBS) \
			-DHAVE_CONFIG_H -I$(PYTHON_INCLUDE)

			$(CC)  -shared pyhsbrd.o pyhsbrd_wrap.o -o _pyhsbrd.so \
			$(LIBS_PATH) $(LIBS)

clean:
	rm -f *.o *~ core *.cxx *.py *.so
