/*
 H2_main.c
Main program for a variational Monte Carlo simulation of a helium atom.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "func.h"
#define PI 3.141592653589

// Main program 
int main(){

	// Declaration of variables and arrays
	int i, j;
	double sum, sum2;
	int N; // Number of interations
	double mean, mean2, var; // <f>, <f^2> and var[f]
	double delta;
	double q;
	int throw_away, norejection; // Number of iterations to throw away in the begining, number of rejections
	double random; // Random number [0,1]
	double alpha;
	double positions[2][3]; // Positions in 3D for 2 particles
	double temp[2][3]; // Temporary array for new positions
	double p, p_temp; // Probabilities
	double distance; 
	double wave_func;
	double energy;
	double distances_nucleus[2];
	
	// Initialize variables
	sum = 0;
	sum2 = 0;
	var = 0;
	delta = 0.65;
	throw_away = 2000;
	norejection = 0;
	alpha = 0.1;
	N = 10000;

	srand(time(NULL));

	// Initialize positions
	for(i = 0; i < 3; i++){
		positions[0][i] = 1.0;
		positions[1][i] = -1.0;
	}

	
	// Get initial distances
	distance = getDistance(positions);
	printf("Dist: %f \n", distance);

	// Get wave function
	wave_func = get_wavefunction(positions, alpha, distance);
	printf("W: %f \n", wave_func);

	// Get wave function
	energy = get_local_e(positions, alpha);
	printf("E: %f \n", energy);

	// Calculate the probability (Not normalized)
	p = pow(wave_func, 2);
	printf("P: %f \n", p);

	// Open a file to print the variable x in
	FILE *m_file;
	m_file = fopen("distances.data","w");

	// Get distances to nucleus
	get_distances_nucleus(positions, distances_nucleus);

	// Save initial distances to nucleus
	fprintf(m_file,"%f \n", distances_nucleus[0]);
	fprintf(m_file,"%f \n", distances_nucleus[1]);
	
	// Calculate the integral
	for(j = 1; j < N; j++){

		// Generate random numbers and get next configuration
		for(i = 0; i < 3; i++){
			random = (double) rand() / (double) RAND_MAX;	
			temp[0][i] = positions[0][i] + delta*(random - 0.5);

			random = (double) rand() / (double) RAND_MAX;	
			temp[1][i] = positions[1][i] + delta*(random - 0.5);
		}

		// Calculate distance between the particles
		distance = getDistance(positions);

		// Get wave function
		wave_func = get_wavefunction(positions, alpha, distance);

		// Calculate the probability
		p_temp = pow(wave_func,2);
		q = p_temp/p;
		
		// If q > 1 all trials will be accepted
		if (q < 1){

			// New random number
			random = (double) rand() / (double) RAND_MAX;

			// Trial, if q >= random, save the temporary positions
			if(q >= random){
				for(i = 0; i < 3; i++){
					positions[0][i] = temp[0][i];
					positions[1][i] = temp[1][i];
				}	
				p = p_temp;
				norejection++;
				}
		}
		else{
			for(i = 0; i < 3; i++){
					positions[0][i] = temp[0][i];
					positions[1][i] = temp[1][i];
			}	
			p = p_temp;
			norejection++;

		}
		// Skip the 'throw_away' first datapoints
		if(j > throw_away){
			//sum += distance * (1-distance) * 2.0 / sin(PI * distance) / PI;
			//sum2 += distance * (1-distance) * 2.0 / sin(PI * distance) / PI * distance * (1-distance) * 2.0 / sin(PI * distance) / PI;
		}

		// Get distances to nucleus
		get_distances_nucleus(positions, distances_nucleus);

		// Save initial distances to nucleus
		fprintf(m_file,"%f \n", distances_nucleus[0]);
		fprintf(m_file,"%f \n", distances_nucleus[1]);

	}

	// Get the means to calculate the variance
	mean = sum/(N-throw_away-1);
	mean2 = sum2/(N-throw_away-1);
	var = (mean2 - mean*mean)/(N-throw_away-1);	
	
	// In the terminal, print how many throw aways
	printf("Nbr of rejections: %i \n", N-norejection);

	// Close the data-file
	fclose(m_file); 

}
