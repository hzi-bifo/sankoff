ALL= build/sankoff

build/sankoff: src/sankoff.cc src/tree.h
	g++ -o build/sankoff src/sankoff.cc -std=c++11 -Wall -lboost_program_options -I/home/sgoliaei/miniconda3/include/ -L/home/sgoliaei/miniconda3/lib/ -ltbb
