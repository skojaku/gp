# Makefile
.PHONY: all

MEXFLAGS := -fopenmp
MEXOFLAGS := -O3
MEXCOMPILER := /usr/bin/gcc-4.9
CC := g++

CFLAGS := -O3 -std=c++11 -fopenmp

# for linux
all: gp_mex.mexa64 gp

gp_mex.mexa64: ./src/*
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' GCC='$(MEXCOMPILER)' src/gp_mex.cpp

gp: ./src/*
	$(CC) $(CFLAGS) -o gp src/gp_cl.cpp 
	

.PHONY: clean
clean:
	$(RM) *.mexa64 
