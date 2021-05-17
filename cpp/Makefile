
test: lib.cpp test.cpp
	g++ -Wall --std=c++11 -O3 -march=native -I. -c -o lib.o lib.cpp 
	g++ -Wall --std=c++11 -O3 -march=native -I. -shared  lib.o -o lib.so
	g++ -Wall --std=c++11 -O3 -march=native -I. -o test test.cpp

exec: test
	./test

clean:
	rm test lib.o lib.so

