Used gcc 7 and newest version of mpicc

COMPILING:
----------

Command to compile the serial version:
'g++ serial.cpp -o example_executable_file -Ofast'

Command to compile the OpenMP version:
'g++ -fopenmp serial.cpp -o example_executable_file -Ofast'


RUNNING:
----------

Command to run the serial version:
'./serial -i example_input_file.txt -o example_output_file.txt

Command to run OpenMP version:
'./openmp -i example_input_file.txt -o example_output_file.txt

