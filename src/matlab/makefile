# Makefile
.PHONY: all

#MEXFLAGS := -fopenmp
MEXOFLAGS := -O3 -fopenmp
#MEXOFLAGS := -Wall -g -fopenmp
MEXCOMPILER := g++
#MEXCOMPILER := /usr/bin/gcc-4.9
CC := g++

CFLAGS := -O3 -std=c++11 -fopenmp

# for linux
all: gp_mex.mexa64 qstest_mex.mexa64

gp_mex.mexa64: ../lib/gp.h ../lib/community-detection-algorithms/* gp_mex.cpp ../lib/quality_functions.h
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' GCC='$(MEXCOMPILER)' gp_mex.cpp

#gp: ../lib/gp.h ./src/community-detection-algorithms/* ./src/gp_cl.cpp ./src/quality_functions.h
#	$(CC) $(CFLAGS) -o gp src/gp_cl.cpp 

qstest_mex.mexa64: ../lib/qstest/* ../lib/community-detection-algorithms/* qstest_mex.cpp ../lib/quality_functions.h
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' GCC='$(MEXCOMPILER)' qstest_mex.cpp
	

.PHONY: clean
clean:
	$(RM) *.mexa64 ./gp  
