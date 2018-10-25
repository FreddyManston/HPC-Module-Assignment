#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "includes/dcdplugin.c"
#include <time.h>

/*
** BRUTE FORCE IMPLEMENTATION OF CLOSEST-PAIR PROBLEM
*/

float getDistance(float *atom1, float *atom2) {
	//printf("%3.2f,	%3.2f,	%3.2f\n", atom1[0], atom1[1], atom1[2]);
	//printf("%3.2f,	%3.2f,	%3.2f\n", atom2[0], atom2[1], atom2[2]);

	float X = atom1[0] - atom2[0];
	float Y = atom1[1] - atom2[1];
	float Z = atom1[2] - atom2[2];

	//printf("%3.2f\n", X*X);
	//printf("%3.2f\n", Y*Y);
	//printf("%3.2f\n", Z*Z);

	return(X*X + Y*Y + Z*Z);
}

/*
	INPUT 1:
		example_pn3_10RU_751frames.dcd
		3
		1-403
		168043-168052
		
	INPUT 2:
		example_pn3_10RU_751frames.dcd
		15
		1-403
		404-168042
*/

int main(int argc, char **argv) {
	/*
	 *	Getting the input from the input file:
	*/
	const char *dcd_file = "example_pn3_10RU_751frames.dcd";
	char *setA = "1-403";
	char *setB = "404-168052";
	const int k = 15;

	int set_A[403];
	int set_B[167639];

	// Creating A
	for (int i = 0; i <= 402; i++) { set_A[i] = i + 1; }
	// Creating B
	for (int i = 0; i <= 167638; i++) { set_B[i] = i + 404; }
	// Size of A and B
	int size_A = sizeof(set_A)/sizeof(set_A[0]);
	int size_B = sizeof(set_B)/sizeof(set_B[0]);

	/*
	 *	Rest of the code:
	*/
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
		/* At this point 
			dcd contains
			dcd->x = Array of X coordinates of all atoms for timestep step
			dcd->y = Array of Y coordinates of all atoms for timestep step
			dcd->z = Array of Z coordinates of all atoms for timestep step
			
			timestep contains
			timestep.coords = Array of packed XYZ coordinates of all atoms for timestep step
							 [X1, Y1, Z1, X2, Y2, Z2, ..., Xn, Yn, Zn] where n = natoms
		  
			Both are overwritten next loop 
		*/

		// Outputs and smallest dist. of atoms:
		int frame[k], atoma[k], atomb[k];
		float k_smallest[k];
		for (int i = 0; i < k; i++) { 
			k_smallest[i] = 10000.00;	// assume largest distance
			frame[i] = step;
			atoma[i] = 0;
			atomb[i] = 0;
		}

		//printf("Timestep %d\n", step);
		//printf("Step:	X	Y	Z\n\n");
		for (int i = 0; i < size_A; ++i) {
			float *current1 = timestep.coords + 3 * (set_A[i]);
			for (int j = 0; j < size_B; ++j) {
				float *current2 = timestep.coords + 3 * (set_B[j]);
				float dist = getDistance(current1, current2);	// No need to cater for dist=0 because A and B are disjoint

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
		//for (int i = 0; i < k; i++) { printf("%d,%d,%d,%f\n", frame[i], atoma[i], atomb[i], sqrt(k_smallest[i])); }
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