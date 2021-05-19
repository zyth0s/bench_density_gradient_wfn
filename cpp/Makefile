test: test.cpp lib.so
	g++ -Wall --std=c++11 -O3 -march=native -fPIC -I. -o test test.cpp

lib.so: lib.o
	g++ -Wall --std=c++11 -O3 -march=native -fPIC -I. -shared  lib.o -o lib.so

lib.o: lib.cpp
	g++ -Wall --std=c++11 -O3 -march=native -fPIC -I. -c -o lib.o lib.cpp

runtest: test
	./test

clean:
	rm test lib.o lib.so

.PHONY: clean runtest
