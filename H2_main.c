/*
 H2_main.c

Main program for a variational Monte Carlo simulation of a helium atom. 
The simulation iterates for a total of "N" iterations. 
It equilibrated the system for "throw_away" number of iterations. 
Results for the equilibrating interations are not written to the outputfile. 
It is also possible to optimize the "alpha"-parameter by rescaling it 
with the function "rescale_alpha". This function is turned on by putting 
the char-variable "rescale_on" to 'y'. For the different tasks in the assignment, 
use the following settings for "*** Variables to change for different tasks ***". 

** TASK 1: ** 
alpha_start = 0.1;
alpha_stop = 0.1;
N = 100000;
throw_away = 0;
rescale_on  = 'n';
rescale_after_iterations = 0;
beta = 0; 

** TASK 2: **
alpha_start = 0.1;
alpha_stop = 0.1;
N = 1000000; 
throw_away = 10000; 
rescale_on  = 'n';
rescale_after_iterations = 0;
beta = 0; 

Plot energy.data and determine how big throw_away should be. 
Set throw_away and make a proper simulation.

** TASK 3: **
alpha_start = 0.05;
alpha_stop = 0.25;
N = 1000000;
throw_away = 10000;
rescale_on  = 'n';
rescale_after_iterations = 0;
beta = 0; 

** TASK 4: **
alpha_start = 0.1;
alpha_stop = 0.1;
N = 1000000;
throw_away = 10000;
rescale_on  = 'y';
rescale_after_iterations = 10000;
beta = 0.8; 

Simulate for different values of beta and plot the data. 
Try for beta = {0.6, 0.7, 0.75, 0.8, 0.9}.

** TASK 5: **
alpha_start = 0.1;
alpha_stop = 0.1;
N = 2000000;
throw_away = 100000;
rescale_on  = 'n';
rescale_after_iterations = 0;
beta = 0; 

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
	int i, j, n, m, k;
	int N; // Number of interations
	double delta; // Correction parameter for generating new configurations
	double q;
	int throw_away, norejection; // Number of iterations to throw away in the begining, number of rejections
	double random; // Random number [0,1]
	double alpha_start, alpha_stop, alpha, new_alpha;
	double positions[2][3]; // Positions in 3D for 2 particles
	double temp[2][3]; // Temporary array for new positions
	double p, p_temp; // Probabilities
	double distance; // Distance between the two electrons
	double wave_func;
	double energy_mean; // Moving average
	double distances_nucleus[2]; // Distance between the electrons and the nucleus
	int iteration; // Iteration number for rescaling alpha
	double alpha_sum;
	char rescale_on;
	int rescale_after_iterations;
	double beta; // Parameter to rescale alpha
	int nbr_simulations;

	// Initialize variables
	delta = 0.967;
	alpha = 0;

	// *** Variables to change for different tasks ***
	alpha_start = 0.10;
	alpha_stop = 0.10;
	N = 200000;
	throw_away = 100000;
	rescale_on  = 'y'; // Rescale alpha for rescale_on = 'y'
	rescale_after_iterations = 1000; // Rescale alpha after rescale_after_iterations
	beta = 0.8; // Shouble be (0.5,1]
	nbr_simulations = 1; // Should be > than 300 if a mean of the mean is wanted, else should be 1

	// Allocate memory for big arrays
	double *energy_l = malloc((N + 1) * sizeof(double));
	double *grad_ln_wave = malloc((N + 1) * sizeof(double));
	double *mean = malloc(nbr_simulations * sizeof(double));

	// Seed for generating random numbers
	srand(time(NULL));

	// Open a file to print the distance in
	FILE *m_file;		
	m_file = fopen("distances.data","w");

	// Open a file to print the energy and the alpha value
	FILE *e_file;
	e_file = fopen("energy.data","w");

	for(k = 0; k < nbr_simulations; k++){

		// Perform the simulation for different alpha values
		for(alpha = alpha_start; alpha <= alpha_stop; alpha = alpha + 0.025){

			// Initiation for each new loop
			norejection = 0;
			energy_mean = 0;
			alpha_sum = 0;
			n = 1;
			m = 1;

			// Print what alpha
			printf("********** ALPHA = %.3f **********\n", alpha);

			// Initialize positions
			for(i = 0; i < 3; i++){
				positions[0][i] = 10.0;
				positions[1][i] = -10.0;
			}

			// Generate random numbers and get small displacements in the initial configuration
			for(i = 0; i < 3; i++){
				random = (double) rand() / (double) RAND_MAX;				
				temp[0][i] = positions[0][i] + delta*(random - 0.5);

				random = (double) rand() / (double) RAND_MAX;	
				temp[1][i] = positions[1][i] + delta*(random - 0.5);
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
			alpha_sum = alpha;

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
					if(rescale_on == 'y' && j%rescale_after_iterations == 0){

						n++;
						new_alpha = rescale_alpha(new_alpha, energy_l, grad_ln_wave, distance, j - throw_away, rescale_after_iterations, beta);		
						alpha_sum += new_alpha;
					}
				}

			}

			// Get statistical inefficiency from the correlation function
			mean[k] = error_corr_func(energy_l, N + 1 - throw_away);

			// Get statistical inefficiency from block averaging
			error_block_average(energy_l, N + 1 - throw_away);

			// In the terminal, print how many rejections
			printf("Tot nbr iteration: %i \nNbr eq iterations: %i \n", N, throw_away);
			printf("Nbr rejections: %i\nMean alpha: %f\n", N-norejection, alpha_sum/n);

			// Print loop finished-line
			printf("***********************************\n");
		}
	
	}

	// If many independent simulation, take a mean and determine the error bar
	if(k > 1){
		printf("For %i simulations: \n", k);
		error_corr_func(mean, k);
	}

	// Close the data-files
	fclose(m_file); 
	fclose(e_file);

	free(energy_l); free(grad_ln_wave); free(mean);
	energy_l = NULL; grad_ln_wave = NULL; mean = NULL;
}

