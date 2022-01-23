all:
	g++ -Wall -flto -O3 -I headers src/*.cpp main.cpp -o bloomfilter
	
lto_off: 
	g++ -Wall -O3 -I headers src/*.cpp main.cpp -o bloomfilter

unoptimized:
	g++ -Wall -I headers src/*.cpp main.cpp -o bloomfilter

clang:
	clang++ -Wall -flto -O3 -I headers src/*.cpp main.cpp -o bloomfilter	

clang_lto_off:
	clang++ -Wall -O3 -I headers src/*.cpp main.cpp -o bloomfilter

clang_unoptimized:
	clang++ -Wall -I headers src/*.cpp main.cpp -o bloomfilter

clean:
	rm -f *.o bloomfilter
