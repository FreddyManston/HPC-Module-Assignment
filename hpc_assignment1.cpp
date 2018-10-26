#include "includes/dcdplugin.c"
#include "includes/molfile_plugin.h"
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
#include <omp.h>
using namespace std;

void print(std::vector<int> const &input)
{
  for (int i = 0; i < input.size(); i++) {
    std::cout << input.at(i) << ' ';
  }
}

class Pair{
public:
	int a,b;
	float dist;

	Pair(const int& indexA, const int& indexB, const float& d){
		a = indexA;
		b = indexB;
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

	/*SEPARATE INPUT BY COMMA*/
	int current, previous = 0;
	current = indices.find(",");
	while(current != std::string::npos){
		separateCommas.push_back(indices.substr(previous, current-previous));
		previous = current+1;
		current = indices.find(",", previous);
	}
	separateCommas.push_back(indices.substr(previous, current-previous));

	/*ADD ALL INDEX VALUE*/
	for (int i = 0; i < separateCommas.size(); i++) {
		if(separateCommas.at(i).find("-") != std::string::npos){
			int setAstart, setAend;
			stringstream helper(separateCommas.at(i).substr(0, separateCommas.at(i).find("-")));helper >> setAstart;
			stringstream helper1 (separateCommas.at(i).substr(separateCommas.at(i).find("-")+1));helper1 >> setAend;
			for (int j = setAstart; j <= setAend; j++) {
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
	/*READ FROM FILE*/
	const auto argI = std::string{"-i"};
	const auto argO = std::string{"-o"};

	std::string inFileName, outFileName;

	for (int i = 0; i < argc; i++) {
		if (argv[i] == argI && i <argc -1){
			inFileName = argv[i+1];
		} else if(argv[i] == argO && i <argc-1){
			outFileName = argv[i+1];
		}
	}
	std::ifstream infile;
	infile.open(inFileName);

	string file_name="";
	string k_most="";string setA="";string setB="";
	int setAstart, setAend, setBstart, setBend, kMin;

	getline(infile, file_name);
	getline(infile, k_most);
	getline(infile, setA);
	getline(infile, setB);

  std::cout << "DCD FILE: " << file_name << std::endl;
  std::cout << "k: " << k_most << std::endl;
  std::cout << "A range: " << setA << std::endl;
  std::cout << "B range: " << setB << std::endl;

	stringstream helper (k_most); helper >> kMin;

	vector<int> setAindices = unroll(setA);
	vector<int> setBindices = unroll(setB);

  std::cout << "Printing atoms in A:\n";
  for (int i = 0; i < setAindices.size(); i++) {
    std::cout << setAindices.at(i) << ' ';
  }std::cout << std::endl;
  std::cout << "Printing atoms in B:\n";
  for (int i = 0; i < setBindices.size(); i++) {
    std::cout << setBindices.at(i) << ' ';
  }std::cout << std::endl;

	ofstream output_file, time_file;
	output_file.open(outFileName); time_file.open("timings.txt");

  cout << "dcd file: " << file_name << std::endl;
  file_name = "example_pn3_10RU_751frames.dcd";
  cout << "dcd file: " << file_name << std::endl;
  
	/*USE DCDPLUGIN TO READ DCD FILE*/
	int natoms = 0;
	void *raw_data = open_dcd_read(file_name.c_str(), "dcd", &natoms);
	if (!raw_data){
    std::cout << "Error reading file: "<< file_name << std::endl;
    return 1;
  }

	dcdhandle *dcd = (dcdhandle *)raw_data;

	cout << "Number of atoms: " << natoms << endl;
	cout << "Number of timesteps: " << dcd -> nsets << endl;

	/*PERFORM K-CLOSEST PAIRS ON EACH TIMESTEP*/
	clock_t begin = clock();
	auto start = chrono::steady_clock::now();

	int globalTimeStep = -1;
	int numTimesteps = dcd->nsets;
	std::priority_queue<Pair, std::vector<Pair>, distCompare> k_dist_timesteps[numTimesteps];

	#pragma omp parallel for
	for (int i = 0; i < numTimesteps; i++) {

		molfile_timestep_t timestep;
		timestep.coords = (float *)malloc(3*sizeof(float)*natoms);

		int localTimeStep;
		#pragma omp critical
		{
			globalTimeStep++;
			localTimeStep = globalTimeStep;
			int rc = read_next_timestep(raw_data, natoms, &timestep);
		}

		std::priority_queue<Pair, std::vector<Pair>, distCompare> k_dist;
		Pair fake(0,0,1000);
		k_dist.push(fake);
		int indexA, indexB;
		float x1,y1,z1,x2,y2,z2,dist;
		for (int j = 0; j < setAindices.size(); j++) {
			indexA = setAindices.at(j);
			x1 = *(timestep.coords + 3*indexA);
			y1 = *(timestep.coords + 3*indexA+1);
			z1 = *(timestep.coords + 3*indexA+2);
			for (int k = 0; k < setBindices.size(); k++) {
				indexB = setBindices.at(k);
				x2 = *(timestep.coords + 3*indexB);
				y2 = *(timestep.coords + 3*indexB+1);
				z2 = *(timestep.coords + 3*indexB+2);
				dist = pow((x1 - x2), 2.0) + pow((y1 - y2), 2.0) + pow((z1 - z2), 2.0);
				if (dist < k_dist.top().dist){
					Pair temp(indexA, indexB, dist);
					k_dist.push(temp);
					if (k_dist.size() > kMin){
						k_dist.pop();
					}
				}
			}
		}
			k_dist_timesteps[localTimeStep] = k_dist;
	}

	for (int i = 0; i < numTimesteps; i++){
		for (int j = 0; j < kMin; j++)
		{
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
