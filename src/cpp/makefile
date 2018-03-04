# Makefile
.PHONY: all

CC := g++

CFLAGS := -O3 -std=c++11 -fopenmp

# for linux
all: gp

gp: ../lib/gp.h ../lib/community-detection-algorithms/* gp_cl.cpp ../lib/quality_functions.h
	$(CC) $(CFLAGS) -o gp gp_cl.cpp 


.PHONY: clean
clean:
	$(RM) ./gp  
