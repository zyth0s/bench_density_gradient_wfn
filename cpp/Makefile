test: test.cpp lib.so
	g++ -Wall --std=c++11 -O3 -march=native -fPIC -funroll-loops -I. -o test test.cpp

lib.so: lib.o
	g++ -Wall --std=c++11 -O3 -march=native -fPIC -funroll-loops -I. -shared  lib.o -o lib.so

lib.o: lib.cpp
	g++ -Wall --std=c++11 -O3 -march=native -fPIC -funroll-loops -I. -c -o lib.o lib.cpp

runtest: test
	./test

clean:
	rm test lib.o lib.so

.PHONY: clean runtest
