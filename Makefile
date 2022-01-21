all:
	g++ -I headers src/*.cpp main.cpp -o bloomfilter

clean:
	rm -f *.o bloomfilter *.core
