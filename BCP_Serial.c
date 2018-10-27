#include "includes/dcdplugin.c"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <string.h>
#include <omp.h>
#include "groups.h"
//using namespace std;
/*
** BRUTE FORCE IMPLEMENTATION OF CLOSEST-PAIR PROBLEM
*/

float distance(float *atom1, float *atom2) {

	float X = atom1[0] - atom2[0];
	float Y = atom1[1] - atom2[1];
	float Z = atom1[2] - atom2[2];

	return(X*X + Y*Y + Z*Z);
}

static int getRange(char *indec){
	//printf("indeces are %s\n",indeces);
	int size=0;
	char l[256];
	strcpy(l,indec);
	char *e;
	int ind;
	e=strchr(l,'-');
	ind=(int) (e-l);
	char left[100];
	char right[100];
	memcpy(left,&l[0],ind);
	memcpy(right,&l[ind]+1,strlen(l));
	//printf("%s and %s \n",left, right );
	//int low=atoi(strtok(NULL,"-"));
	//int hi=atoi(strtok(NULL,","));
	//printf("found a range from %s to %s \n",left,right);
	//printf("range is %i\n", atoi(right)-atoi(left));
	return (atoi(right)-atoi(left)+1);//(upper-lower);

}

//returns atom indeces in an array
int * getAts(char *ind,int size){
	int *atoms = malloc(sizeof(int)*size);
	int counter = 0;
	char line[256];
	ind[strcspn(ind, "\r\n")] = 0;
	strcpy(line, ind);
	char *token = strtok(line,",");
	while (token != NULL) {
		if (strchr(token,'-')!=NULL) {
			//printf("dash foune %s\n",token);
			int size = 0;
			char l[256];
			strcpy(l,token);
			char *e;
			int ind;
			e = strchr(l,'-');
			ind = (int) (e-l);
			char left[100];
			char right[100];
			memcpy(left,&l[0],ind);
			memcpy(right,&l[ind]+1,strlen(l));
			for (int s = atoi(left); s<=atoi(right); s++) {
				atoms[counter] = s;
				counter = counter + 1;	
			}
		}
		else {
			atoms[counter] = atoi(token);
			counter = counter + 1;
		}
		token = strtok(NULL,",");
		//printf("line is %s token is %s \n",line,token);
	}
	return atoms;
}

struct group * getList(char *indeces) {
	int groupSize = 0; // we predetermine the size of the array needed to hold the atom indexes
	char line[256];
	indeces[strcspn(indeces, "\r\n")] = 0;
	strcpy(line, indeces);
	char *token = strtok(line,",");
	while (token != NULL) {
		if (strchr(token,'-') != NULL) {
			printf("dash found %s\n",token);
			groupSize = groupSize + getRange(token);
		}
		else {
			printf("beat me up %s\n",token );
			groupSize = groupSize+1;
		}
		token = strtok(NULL,",");
		printf("line is %s token is %s \n",line,token);
	}
	int *atoms = malloc(sizeof(int)*groupSize);
	atoms = getAts(indeces,groupSize);
	struct group* g = (struct group*)malloc(sizeof(struct group));
	g->groupSize = groupSize;
	g->groupAtoms = atoms;
	//free(atoms);
	return g;
}

int main(int argc, char *argv[]) {
	///////////////// reading from file //////////////////
	int getopt(int argc, char * const argv[], const char *optstring);
	extern char *optarg;
	extern int optind, opterr, optopt;
	int ich;
	char *fileName;
	while ((ich = getopt (argc, argv, "ai:c")) != EOF) {
		switch (ich) {
			case 'a':	// Flags/Code when -a is specified
				break;
			case 'i': 
				fileName=optarg;
				printf("found file name %s\n",fileName);
				// Flags/Code when -b is specified
				// The argument passed in with b is specified
				// by optarg
				break;
			case 'c': // Flags/Code when -c is specified
				break;
			default: // Code when there are no parameters
				printf("default file name\n");
				fileName="input";
				break;
		}
	}
	
	char dcd_file[256];
	char tmp[256];
	char a[256];
	char b[256];
	int k;
  	//fileName="example_input_file2.txt";
  	//printf("files name is %s\n",fileName);
	if(fileName!=NULL){
		FILE *file = fopen ( fileName, "r" );
		//printf("got here\n");
		if ( file != NULL ){
			if(fgets(dcd_file, 100, file) && fgets(tmp, 100, file) && fgets(a, 100, file) && fgets(b, 100, file)){
				dcd_file[strcspn(dcd_file,"\r\n")]=0;
				tmp[strcspn(tmp,"\r\n")]=0;
				a[strcspn(a,"\r\n")]=0;
				b[strcspn(b,"\r\n")]=0;
				k=atoi(tmp);
				//printf("dcd is %s %s %s %s\n",dcd_file,tmp,a,b);	
				fclose ( file );
			}
			else{
				printf("input error in text file");
				return 0;	
			}	
   		}
   		else{
	  		printf("file named %s was not found\n",fileName);
	  		return 0;
   		}
		
	}
	else{
		printf("file named %s is empty?\n",fileName);
		return 0;
	}

	//now we call a function to return group structures,
	//parsing the indeces from the textfile
	struct group *g1 =getList(a);
	//printf("size is  %d\n", g1->groupSize);
	struct group *g2 =getList(b);
	printf("Here now\n");
	printf("Here's all the input: %s, %s, %s, %d\n", dcd_file, a, b, k);


////////////////////////////////////////////////////////////////////////////////////////

	int set_A[403];
	int set_B[167639];

	// Creating A
	for (int i = 0; i <= 402; i++) { set_A[i] = i + 1; }
	// Creating B
	for (int i = 0; i <= 167638; i++) { set_B[i] = i + 404; }
	// Size of A and B
	int size_A = sizeof(set_A)/sizeof(set_A[0]);
	int size_B = sizeof(set_B)/sizeof(set_B[0]);

	//// Rest of the code
	int natoms = 0;
	void *raw_data = open_dcd_read(dcd_file, "dcd", &natoms);
	if (!raw_data) { return 1; }

	dcdhandle *dcd = (dcdhandle *) raw_data;

	molfile_timestep_t timestep;
	timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);
	int read_failed = 0;
	int rc;

	clock_t begin = clock();	// Start timer
	for (int step = 0; step < dcd->nsets; ++step) {
		rc = read_next_timestep(raw_data, natoms, &timestep);
		if (rc) {
			read_failed = 1;
			break;
		}

		// Outputs and smallest dist. of atoms:
		int frame[k], atoma[k], atomb[k];
		float k_smallest[k];
		for (int i = 0; i < k; i++) { 
			k_smallest[i] = 10000.00;	// assume largest distance
			frame[i] = step;
			atoma[i] = 0;
			atomb[i] = 0;
		}

		for (int i = 0; i < size_A; ++i) {
			float *current1 = timestep.coords + 3 * (set_A[i]);
			for (int j = 0; j < size_B; ++j) {
				float *current2 = timestep.coords + 3 * (set_B[j]);
				float dist = distance(current1, current2);	// No need to cater for dist=0 because A and B are disjoint

				// Checking for and inserting a new smallest distance
				if (dist < k_smallest[k -1]) {
					k_smallest[k - 1] = dist;
					atoma[k - 1] = set_A[i];
					atomb[k - 1] = set_B[j];

					int index = k - 2;
					while (index >= 0 && dist < k_smallest[index]) {
						int temp_atoma = atoma[index];
						int temp_atomb = atomb[index];
						float temp_smallest = k_smallest[index];

						atoma[index] = atoma[index + 1];
						atomb[index] = atomb[index + 1];
						k_smallest[index] = dist;

						k_smallest[index + 1] = temp_smallest;
						atoma[index + 1] = temp_atoma;
						atomb[index + 1] = temp_atomb;

						index --;
					}
				}
			}
		}
		for (int i = 0; i < k; i++) { printf("%d,%d,%d,%f\n", frame[i], atoma[i], atomb[i], sqrt(k_smallest[i])); }
		if (step >= 751) break;
	}

	clock_t end = clock();	// End timer
	float time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Total time taken: %f\n", time_spent);

	free(timestep.coords);
	close_file_read(raw_data);
	if (read_failed) {
		return 2;
	}

	printf("DONE\n");

	return 0;
}
