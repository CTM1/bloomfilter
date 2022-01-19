all:
	g++ -o bloomfilter main.cpp

clean:
	rm -f *.o bloomfilter *.core
