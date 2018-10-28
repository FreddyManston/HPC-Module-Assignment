Used gcc 7 and newest version of mpicc


THIS REPO CONTAINS THE FOLLOWING:
---------------------------------

data/			- the inputs used in connection with the dcd file
				- example outputs

includes/ 		- neccessary files for reading the dcd file

BCP_Serial.c 	- The initial. completed, attempt at the problem, in C

makefile		- run 'make' to compile all relevant code

mpi.cc 			- an attempt at the MPI implementation

serial.cc 		- both the serial AND OpenMP versions, compiled with or without -fopenmp flag


COMPILING:
----------

Command to compile the serial version:
'g++ serial.cpp -o serial.o -Ofast'

Command to compile the OpenMP version:
'g++ -fopenmp serial.cpp -o openmp.o -Ofast'


RUNNING:
----------

Command to run the serial version:
'./serial.o -i example_input_file.txt -o example_output_file.txt

Command to run OpenMP version:
'./openmp.o -i example_input_file.txt -o example_output_file.txt

