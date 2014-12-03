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
	int throw_away, rejections; // Number of iterations to throw away in the begining, number of rejections
	double random[2][3]; // Random number [0,1]
	double alpha;
	double positions[2][3]; // Positions in 3D for 2 particles
	double p[2][3]; // Probabilities

	// Initialize variables
	sum = 0;
	sum2 = 0;
	var = 0;
	delta = 0.45;
	throw_away = 2000;
	rejections = 0;
	alpha = 0.1;
	N = 10000;

	// Initialize positions
	for(i = 0; i < 3; i++){
		positions[0][i] = 1.0;
		positions[1][i] = -1.0;

		p[0][i] = 0;
		p[1][i] = 0;
	}

	// Open a file to print the variable x in
	FILE *m_file;
	m_file = fopen("positions.data","w");

	// Save initial positions
	fprintf(m_file,"%F \t %F \t %F  \t %F \t %F \t %F \n", positions[0][0], positions[0][1], positions[0][2], positions[1][0], positions[1][1], positions[1][2]);

	// Calculate the integral
	for(j = 1; j < N; j++){

		// Generate random number and get next configuration
		for(){
			random = (double) rand() / (double) RAND_MAX;	
			positions[][j] = x[j-1] + delta*(random - 0.5);
		}

		// Calculate the probability
		p[j] = 0.0;
		p[j] = sin(PI * x[j]);
		q = p[j]/p[j-1];

		// New random number
		random = (double) rand() / (double) RAND_MAX;

		// Trial
		if(q < r){
			x[j] = x[j-1];
			p[j] = p[j-1];
			rejections++;
		}
		
		// Skip the 'throw_away' first datapoints
		if(j > throw_away){
			n++;
			sum += x[j] * (1-x[j]) * 2.0 / sin(PI * x[j]) / PI;
			sum2 += x[j] * (1-x[j]) * 2.0 / sin(PI * x[j]) / PI * x[j] * (1-x[j]) * 2.0 / sin(PI * x[j]) / PI;
		}

	}

	// Get the means to calculate the variance
	mean = sum/(N-throw_away-1);
	mean2 = sum2/(N-throw_away-1);
	var = (mean2 - mean*mean)/(N-throw_away-1);	
	
	// In the terminal, print how many throw aways
	printf("Nbr of rejections: %i \n", rejections);

	// Print the result in the terminal
	printf("For N = %i \t Integral = %.8F ± %.8F \n", N-throw_away, mean, sqrt(var));

	// Print x to distribution.data
	for(j = 0; j < N; j++){
		fprintf(m_file,"%F \n", x[j]);
	}

	// Close the data-file
	fclose(m_file); 

}
