export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib/
mkdir tests/out
build/sankoff --cost tests/data/cost-1.txt --aln tests/data/aln1.txt --tree tests/data/tree1.txt
build/sankoff --cost tests/data/cost-1.txt --aln tests/data/aln1.txt --tree tests/data/tree1.txt --out-as tests/out/tree-1-no-dep-asr.txt --out-tree tests/out/tree-1-no-dep-tree.txt --ilabel inode 
