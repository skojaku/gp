# Makefile
.PHONY: all

MEXFLAGS := -fopenmp
MEXOFLAGS := 
#MEXOFLAGS := -O3
CC := g++

#CFLAGS := -O3 -std=c++11 # use this option in case openmp is not working 
CFLAGS := -O3 -std=c++11 -fopenmp

# for linux
all: gp_mex.mexa64 gp

gp_mex.mexa64: ./src/*
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' src/gp_mex.cpp

gp: ./src/*
	$(CC) $(CFLAGS) -o gp src/gp_cl.cpp 
	

.PHONY: clean
clean:
	$(RM) *.mexa64 
