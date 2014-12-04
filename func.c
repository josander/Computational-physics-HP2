/*
func.c
Contains functions for homeproblem 2/b

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589

// Function that returns a distance between two particles
double getDistance(double positions [][3]){
	
	double distance = 0;
	int i;

	// Get the differences squared in each dimension
	for(i = 0; i < 3; i++){
		distance += pow(positions[0][i] - positions[1][i], 2);
	}
	
	return (sqrt(distance));	


}

double get_wavefunction(double *positions, double alpha, double distance){

	int i;
	double r1;
	double r2;
	double wave_func;

	for(i = 0; i < 3; i++){
		r1 += positions[0][i];
		r2 += positions[1][i];
	}

	r1 = sqrt(r1);
	r2 = sqrt(r2);
	
	wave_func = exp(-2 * r1) * exp(-2 * r2) * exp(distance/(2 * (1 + alpha * distance)));

	return wave_func;
	
}

