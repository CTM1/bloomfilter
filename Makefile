all:
	g++ -Wall -O3 -I headers src/*.cpp main.cpp -o bloomfilter

unoptimized: 
	g++ -Wall -I headers src/*.cpp main.cpp -o bloomfilter

clean:
	rm -f *.o bloomfilter *.core
