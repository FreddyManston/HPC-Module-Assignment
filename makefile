default: all

all: serial openmp

serial: serial.cpp
	g++ serial.cpp -o serial.o -Ofast

openmp: serial.cpp
	g++ -fopenmp serial.cpp -o openmp.o -Ofast

clean veryclean:
	rm -f serial.o openmp.o