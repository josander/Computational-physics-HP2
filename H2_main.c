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
	int i, j, n;
	int N; // Number of interations
	double mean, mean2, var; // <f>, <f^2> and var[f]
	double delta;
	double q;
	int throw_away, norejection; // Number of iterations to throw away in the begining, number of rejections
	double random; // Random number [0,1]
	double alpha_start, alpha_stop, alpha, new_alpha;
	double positions[2][3]; // Positions in 3D for 2 particles
	double temp[2][3]; // Temporary array for new positions
	double p, p_temp; // Probabilities
	double distance; 
	double wave_func;
	double energy_mean;
	double distances_nucleus[2];
	int iteration; // Iteration number fot rescaling alpha
	
	// Initialize variables
	var = 0;
	delta = 0.967;
	alpha = 0;
	alpha_start = 0.15;
	alpha_stop = 0.15;
	N = 500000;
	throw_away = 50000;
	n = 0;


	// Allocate memory for big arrays
	double *energy_l = malloc((N + 1) * sizeof(double));
	double *grad_ln_wave = malloc((N + 1) * sizeof(double));

	// Seed for generating random numbers
	srand(time(NULL));

	// Open a file to print the distance in
	FILE *m_file;		
	m_file = fopen("distances.data","w");

	// Open a file to print the energy and the alpha value
	FILE *e_file;
	e_file = fopen("energy.data","w");

	// Perform the simulation for different alpha values
	for(alpha = alpha_start; alpha <= alpha_stop; alpha = alpha + 0.025){

		// Initiation for each new loop
		norejection = 0;
		energy_mean = 0;

		// Print what alpha
		printf("********** ALPHA = %.3f **********\n", alpha);

		// Initialize positions
		for(i = 0; i < 3; i++){
			positions[0][i] = 1.0;
			positions[1][i] = -1.0;
		}



		// Get initial distances
		distance = getDistance(positions);

		// Get wave function
		wave_func = get_wavefunction(positions, alpha, distance);

		// Calculate the probability (Not normalized)
		p = pow(wave_func, 2);

		// Get distances to nucleus
		get_distances_nucleus(positions, distances_nucleus);

		// Save initial distances to nucleus + energy
		fprintf(m_file,"%f \n", distances_nucleus[0]);
		fprintf(m_file,"%f \n", distances_nucleus[1]);

		// Initiate new_alpha
		new_alpha = alpha;

		// Main for-loop
		for(j = 1; j < N + 1; j++){

			// Generate random numbers and get next configuration
			for(i = 0; i < 3; i++){
				random = (double) rand() / (double) RAND_MAX;	
				temp[0][i] = positions[0][i] + delta*(random - 0.5);

				random = (double) rand() / (double) RAND_MAX;	
				temp[1][i] = positions[1][i] + delta*(random - 0.5);
			}

			// Calculate distance between the particles
			distance = getDistance(temp);

			// Get wave function
			wave_func = get_wavefunction(temp, new_alpha, distance);

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
			}else{
				for(i = 0; i < 3; i++){
						positions[0][i] = temp[0][i];
						positions[1][i] = temp[1][i];
				}	
				p = p_temp;
				norejection++;

			}

			// Calculate distance between the particles
			distance = getDistance(positions);

			// Get the gradient of ln(wavefunction) with respect to alpha
			grad_ln_wave[j] = get_grad_ln_wave(distance, new_alpha);

			// Skip the 'throw_away' first datapoints
			if(j > throw_away){
			
				n++;
				// Get energies for the current configuration
				energy_l[j - throw_away - 1] = get_local_e(positions, new_alpha);
				energy_mean += energy_l[j - throw_away - 1];

				// Get the gradient of ln(wavefunction) with respect to alpha
				grad_ln_wave[j - throw_away - 1] = get_grad_ln_wave(distance, new_alpha);

				// Get distances to nucleus
				get_distances_nucleus(positions, distances_nucleus);

				// Save distances to nucleus
				fprintf(m_file,"%f \n", distances_nucleus[0]);
				fprintf(m_file,"%f \n", distances_nucleus[1]);

				// Save current energies
				fprintf(e_file,"%F \t %F \t %F \n", energy_l[j - throw_away - 1], energy_mean/(j - throw_away), new_alpha);

				// Rescale alpha
				new_alpha = rescale_alpha(new_alpha, energy_l, grad_ln_wave, distance, j - throw_away);

			}

			// For each 5000nd iteration, print 
			if(j%50000 == 0){
				printf("%i out of %i steps\n", j, N);
				printf("Alpha: %f \t terms: %i \t index: %i \n", new_alpha, n, j - throw_away-1);
			}
		
		}

		// Get statistical inefficiency from the correlation function
		error_corr_func(energy_l, N + 1 - throw_away);

		// Get statistical inefficiency from block averaging
		error_block_average(energy_l, N + 1 - throw_away);

		// In the terminal, print how many rejections
		printf("Tot nbr iteration: %i \nNbr eq iterations: %i \n", N, throw_away);
		printf("Nbr rejections: %i \nLast alpha: %f\n", N-norejection, new_alpha);

		// Print loop finished-line
		printf("***********************************\n");
	}

	// Close the data-files
	fclose(m_file); 
	fclose(e_file);

	free(energy_l); free(grad_ln_wave);
	energy_l = NULL; grad_ln_wave = NULL;
}

