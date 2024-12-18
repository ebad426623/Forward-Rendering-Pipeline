all: clean rasterizer
clean: 
	rm -f rasterizer
rasterizer: 
	g++ *.cpp -o rasterizer -std=c++11 -O3
