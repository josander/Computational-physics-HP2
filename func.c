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


// Function that calculates the local energy given position and alfa-parameter
double get_local_e(double position[][3], double alfa){

	int i;
	double distance = getDistance(position);
	double alfaDist = 1 + alfa*distance;
	double distInv = 1/distance;
	double alfaDistInv = 1/alfaDist;
	double vector = 0;	
	double magnitude[2] = {0,0};	

	for(i = 0; i < 3; i++){
		magnitude[0] += position[0][i]*position[0][i];
		magnitude[1] += position[1][i]*position[1][i];
	}

	magnitude[0] = 1/sqrt(magnitude[0]);
	magnitude[1] = 1/sqrt(magnitude[1]);

	for(i = 0; i < 3; i++){
		vector += (position[0][i]*magnitude[0] - position[1][i]*magnitude[1])*(position[0][i] - position[1][i]);
	}

	double E = -4 + vector*distInv*pow(alfaDistInv,2) - distInv*pow(alfaDistInv,3) -0.25*pow(alfaDistInv,4) + distInv;
	return(E);  		 

}


// Function that calculates the local wave function
double get_wavefunction(double positions[][3], double alpha, double distance){

	int i;
	double r1;
	double r2;
	double wave_func;

	for(i = 0; i < 3; i++){
		r1 += pow(positions[0][i],2);
		r2 += pow(positions[1][i],2);
	}

	r1 = sqrt(r1);
	r2 = sqrt(r2);
	
	wave_func = exp(-2 * r1) * exp(-2 * r2) * exp(distance/(2 * (1 + alpha * distance)));

	return wave_func;

}

// Get the distances to the nucleus
void get_distances_nucleus(double positions[][3], double distances_nucleus[]){
		
	int i;
	distances_nucleus[0] = 0;
	distances_nucleus[1] = 0;

	for(i = 0; i < 3; i++){
		distances_nucleus[0] += pow(positions[0][i],2);
		distances_nucleus[1] += pow(positions[1][i],2);
	}

	distances_nucleus[0] = sqrt(distances_nucleus[0]);
	distances_nucleus[1] = sqrt(distances_nucleus[1]);
}

