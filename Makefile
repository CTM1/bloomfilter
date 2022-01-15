all:
	g++ -o bloomfilter bloomfilter.cpp

clean:
	rm -f *.o bloomfilter *.core
