include ~/example.mk

CC=mpic++

LDIR =

CFLAGS += -fPIC

OBJ_ENZ = openfpm_core.o

all: openfpm_core

debug: CFLAGS += -ggdb
debug: all

openfpm_core:
			swig -python -c++ openfpm_core.i
			$(CC)  -O3 -g -c --std=c++11  openfpm_core.cpp \
			$< $(INCLUDE_PATH) $(CFLAGS) $(LIBS_PATH) $(LIBS) \
			-DHAVE_CONFIG_H -I -I$(PYTHON_INCLUDE)

			$(CC)  -O3 -g -c --std=c++11 openfpm_core_wrap.cxx \
			$< $(INCLUDE_PATH) $(CFLAGS) $(LIBS_PATH) $(LIBS) \
			-DHAVE_CONFIG_H -I$(PYTHON_INCLUDE)

			$(CC)  -shared openfpm_core.o openfpm_core_wrap.o -o _openfpm_core.so \
			$(LIBS_PATH) $(LIBS)

clean:
	rm -f *.o *~ core *.cxx openfpm_core.py *.so
