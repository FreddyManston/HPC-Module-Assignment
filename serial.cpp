#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <ctime>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <vector>
#include <string.h>
#include <chrono>
#include <functional>
#include "includes/dcdplugin.c"
#include <omp.h>
using namespace std;

class Pair{
public:
	int a,b;
	float dist;

	Pair(const int& atom_A, const int& atom_B, const float& d){
		a = atom_A;
		b = atom_B;
		dist = d;
	}
	string printPair(){
		std::string returnA = std::to_string(a);
		std::string returnB = std::to_string(b);
		std::string returnDist = std::to_string(sqrt(dist));

		return returnA +", "+ returnB +", "+ returnDist;
	}
};

class distCompare{
public:
	bool operator() (const Pair& lhs, const Pair& rhs) const{
			return (lhs.dist < rhs.dist);
		}
};

vector<int> unroll(string indices){
	vector<string> separateCommas;
	vector<int> allIndices;

	/// SEPARATE INPUT BY COMMA
	int current, previous = 0;
	current = indices.find(",");
	while(current != std::string::npos){
		separateCommas.push_back(indices.substr(previous, current-previous));
		previous = current+1;
		current = indices.find(",", previous);
	}
	separateCommas.push_back(indices.substr(previous, current-previous));

	/// ADD ALL INDEX VALUE
	for (int i = 0; i < separateCommas.size(); i++) {
		if(separateCommas.at(i).find("-") != std::string::npos){
			int str_set_Astart, str_set_Aend;
			stringstream helper(separateCommas.at(i).substr(0, separateCommas.at(i).find("-")));helper >> str_set_Astart;
			stringstream helper1 (separateCommas.at(i).substr(separateCommas.at(i).find("-")+1));helper1 >> str_set_Aend;
			for (int j = str_set_Astart; j <= str_set_Aend; j++) {
				allIndices.push_back(j);
			}
		} else {
			int index;
			stringstream helper(separateCommas.at(i)); helper >> index;
			allIndices.push_back(index);
		}
	}
	return allIndices;
}

int main(int argc, char const *argv[]){
///////////////// READING FROM FILE //////////////////
	const auto argI = std::string{"-i"};
	const auto argO = std::string{"-o"};

	std::string input_file, out_file_name;

	for (int i = 0; i < argc; i++) {
		if (argv[i] == argI && i <argc -1){
			input_file = argv[i+1];
		} else if(argv[i] == argO && i <argc-1){
			out_file_name = argv[i+1];
		}
	}
	std::ifstream reader;
	reader.open(input_file);	// Open input file

	string file_name = "", str_k = "", str_set_A = "", str_set_B = "";	// for storing the input from the input file
	/// Getting the input from the input file:
	getline(reader, file_name);	// name of the dcd file
	getline(reader, str_k);	// k value in string
	getline(reader, str_set_A);	// atom set A
	getline(reader, str_set_B);	// atom set B	

	reader.close();		// Close input file

	/// Printing out the input that was taken from the input file:
	std::cout << "DCD FILE: " << file_name << std::endl;
	std::cout << "k: " << str_k << std::endl;
	std::cout << "A range: " << str_set_A << std::endl;
	std::cout << "B range: " << str_set_B << std::endl;

///////////////// PROCESSING THE INPUT //////////////////

	int k;	// the k smallest distances are to be found
	stringstream helper (str_k); helper >> k;	// parse string to int

	// vector (analogous to list) containing the atoms in set_A and set_B, respectively
	vector<int> set_A = unroll(str_set_A);
	vector<int> set_B = unroll(str_set_B);

	// Printing out atoms in set_A and set_B, respectively
	/*std::cout << "Printing atoms in A:\n";
	for (int i = 0; i < set_A.size(); i++) {
		std::cout << set_A.at(i) << ' ';
	} std::cout << std::endl;
	std::cout << "Printing atoms in B:\n";
	for (int i = 0; i < set_B.size(); i++) {
		std::cout << set_B.at(i) << ' ';
	} std::cout << std::endl;*/

	// Opening the output files
	// N.B. must be done before parallelisation to prevent many processes opening the same file
	ofstream output_file, time_file;
	output_file.open(out_file_name); time_file.open("timings.txt");

	cout << "dcd file: " << file_name << std::endl;

	// remove the weird char at the end of the dcd file name line
	// it's not a newline char. not sure what it is. but i dont like it
	if (!file_name.empty() && file_name[file_name.length()-1] != 'd') { 
		file_name.erase(file_name.length()-1); 
	}

///////////////// USE DCDPLUGIN TO READ DCD FILE //////////////////
	int natoms = 0;
	void *raw_data = open_dcd_read(file_name.c_str(), "dcd", &natoms);
	if (!raw_data){
		std::cout << "Error reading file: "<< file_name << std::endl;
		return 1;
	}

	dcdhandle *dcd = (dcdhandle *)raw_data;

	cout << "Number of atoms: " << natoms << endl;
	cout << "Number of timesteps: " << dcd -> nsets << endl;

///////////////// START ACTUAL PROGRAM /////////////////
	clock_t begin = clock();
	auto start = chrono::steady_clock::now();

	int global_step = -1;
	int numTimesteps = dcd->nsets;
	std::priority_queue<Pair, std::vector<Pair>, distCompare> k_dist_timesteps[numTimesteps];	// array of priority queues to store each time step
	
	#pragma omp parallel for
	for (int t = 0; t < numTimesteps; t++) {
		// Inside so that each process has its own instance, i.e. prevents data race
		molfile_timestep_t timestep;
		timestep.coords = (float *)malloc(3*sizeof(float)*natoms);

		int local_step;
		#pragma omp critical
		{
			global_step++;
			local_step = global_step;
			int rc = read_next_timestep(raw_data, natoms, &timestep);
		}

		std::priority_queue<Pair, std::vector<Pair>, distCompare> k_dist;
		Pair fake(0,0,1000);	// assume largest distance
		k_dist.push(fake);
		int atom_A, atom_B;
		float x1,y1,z1,x2,y2,z2,dist;
		for (int i = 0; i < set_A.size(); i++) {
			// Atom # for A set:
			atom_A = set_A.at(i);
			// x, y and z coords for atom A:
			x1 = *(timestep.coords + 3*atom_A);
			y1 = *(timestep.coords + 3*atom_A+1);
			z1 = *(timestep.coords + 3*atom_A+2);
			for (int j = 0; j < set_B.size(); j++) {
				// Atom # for B set:
				atom_B = set_B.at(j);
				// x, y and z coords for atom B:
				x2 = *(timestep.coords + 3*atom_B);
				y2 = *(timestep.coords + 3*atom_B+1);
				z2 = *(timestep.coords + 3*atom_B+2);
				// Getting distance between the two points:
				dist = pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0) + pow((z1 - z2), 2.0);
				if (dist < k_dist.top().dist){
					Pair temp(atom_A, atom_B, dist);
					k_dist.push(temp);
					if (k_dist.size() > k){
						k_dist.pop();
					}
				}
			}
		}
			// free up space, since each process creates its own instance
			// and hex cluster complains
			free(timestep.coords);
			k_dist_timesteps[local_step] = k_dist;
	}

////////// OUTPUTTING EVERYTHING TO FILE //////////
	/// time steps printed in descending order, bc priority queue
	for (int i = 0; i < numTimesteps; i++) {
		for (int j = 0; j < k; j++) {
			Pair top = k_dist_timesteps[i].top();
			k_dist_timesteps[i].pop();
			output_file << i << ", ";
			output_file << top.printPair();
			output_file << "\n";
		}
	}

	clock_t end = clock();
	auto finish = chrono::steady_clock::now();

	double elapsed_secs = double(end-begin) / CLOCKS_PER_SEC;

	// time_file << "real time: " << (chrono::duration_cast<chrono::milliseconds>(finish-start).count())/(double)1000 << "\n";
	// time_file << "cpu time: " << elapsed_secs << "\n";

	cout << "real time: " << (chrono::duration_cast<chrono::milliseconds>(finish-start).count())/(double)1000 << endl;
	cout << "cpu time: " << elapsed_secs << endl;

	close_file_read(raw_data);
	output_file.close();
	time_file.close();

	return 0;
}
