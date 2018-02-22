# Makefile
.PHONY: all

MEXFLAGS := -fopenmp
MEXOFLAGS := -O3
MEXCOMPILER := /usr/bin/gcc-4.9
CC := g++

CFLAGS := -O3 -std=c++11 -fopenmp

# for linux
all: gp_mex.mexa64 gp qstest_mex.mexa64

gp_mex.mexa64: ./src/gp.h ./src/community-detection-algorithms/* ./src/gp_mex.cpp ./src/quality_functions.h
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' GCC='$(MEXCOMPILER)' src/gp_mex.cpp

gp: ./src/gp.h ./src/community-detection-algorithms/* ./src/gp_cl.cpp ./src/quality_functions.h
	$(CC) $(CFLAGS) -o gp src/gp_cl.cpp 

qstest_mex.mexa64: ./src/qstest/* ./src/qstest_mex.cpp ./src/community-detection-algorithms/* ./src/quality_functions.h
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' GCC='$(MEXCOMPILER)' src/qstest_mex.cpp
	

.PHONY: clean
clean:
	$(RM) *.mexa64 ./gp qstest_mex.mexa64  
