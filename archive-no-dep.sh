rm -rf archive/$1.tgz
tar czf archive/$1.tgz src/sankoff.h src/tree.h src/sankoff-no-dep.cc src/GetPot* tests/ data/ makefile
