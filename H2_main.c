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
	double p[2][3]; // Probabilities
	double distance; 

	// Initialize variables
	sum = 0;
	sum2 = 0;
	var = 0;
	delta = 0.45;
	throw_away = 2000;
	norejection = 0;
	alpha = 0.1;
	N = 10000;

	// Initialize positions
	for(i = 0; i < 3; i++){
		positions[0][i] = 1.0;
		positions[1][i] = -1.0;

		p[0][i] = 0;
		p[1][i] = 0;
	}

	// Get initial distances
	distance = 1;//getDistance(positions);

	// Open a file to print the variable x in
	FILE *m_file;
	m_file = fopen("distance.data","w");

	// Save initial positions
	fprintf(m_file,"%F \n", distance);

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
		distance = 1;//getDistance(positions);

		// Calculate the probability
		p[0][j] = sin(PI * distance);
		q = p[0][j]/p[0][j-1];

		// New random number
		random = (double) rand() / (double) RAND_MAX;

		// Trial, if q >= random, save the temporary positions
		if(q >= random){
			for(i = 0; i < 3; i++){
				positions[0][i] = temp[0][i];
				positions[1][i] = temp[1][i];
			}

			p[0][j] = p[0][j-1];
			norejection++;
		}
		
		// Skip the 'throw_away' first datapoints
		if(j > throw_away){
			sum += distance * (1-distance) * 2.0 / sin(PI * distance) / PI;
			sum2 += distance * (1-distance) * 2.0 / sin(PI * distance) / PI * distance * (1-distance) * 2.0 / sin(PI * distance) / PI;
		}

	}

	// Get the means to calculate the variance
	mean = sum/(N-throw_away-1);
	mean2 = sum2/(N-throw_away-1);
	var = (mean2 - mean*mean)/(N-throw_away-1);	
	
	// In the terminal, print how many throw aways
	printf("Nbr of rejections: %i \n", N-norejection);

	// Print the result in the terminal
	printf("For N = %i \t Integral = %.8F ± %.8F \n", N-throw_away, mean, sqrt(var));

	// Print x to distribution.data
	for(j = 0; j < N; j++){
		fprintf(m_file,"%F \n", distance);
	}

	// Close the data-file
	fclose(m_file); 

}
