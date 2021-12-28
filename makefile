ALL= build/sankoff

build/sankoff: src/sankoff.cc src/tree.h
	${CXX} -o build/sankoff src/sankoff.cc -std=c++11 -Wall -lboost_program_options -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -ltbb
