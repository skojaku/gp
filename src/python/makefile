# Makefile
.PHONY: gp

gp:gp.cpp ../lib/*
	clang++ -O3 -Wall -shared -std=c++14 -fPIC `python3 -m pybind11 --includes` gp.cpp -o gp`python3-config --extension-suffix` 
	
.PHONY: clean
clean:
	rm -rf gp 
