export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib/
build/sankoff --cost tests/data/cost-1.txt --aln tests/data/aln1.txt --tree tests/data/tree1.txt
