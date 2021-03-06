/*
func.c
Contains functions for homeproblem 2/b

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589

// Function that returns a distance between the two particles
double getDistance(double positions [][3]){
	
	double distance = 0;
	int i;

	// Get the differences squared in each dimension
	for(i = 0; i < 3; i++){
		distance += pow(positions[0][i] - positions[1][i], 2);
	}
	
	return (sqrt(distance));	


}


// Function that calculates the local energy given position and alpha-parameter
double get_local_e(double position[][3], double alpha){

	int i;
	double distance = getDistance(position);
	double alphaDist = 1 + alpha*distance;
	double distInv = 1/distance;
	double alphaDistInv = 1/alphaDist;
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

	double E = -4 + vector*distInv*pow(alphaDistInv,2) - distInv*pow(alphaDistInv,3) -0.25*pow(alphaDistInv,4) + distInv;
	return(E);
}  		 


// Function that calculates the local wave function
double get_wavefunction(double positions[][3], double alpha, double distance){

	int i;
	double r1 = 0.0;
	double r2 = 0.0;
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
	distances_nucleus[0] = 0.0;
	distances_nucleus[1] = 0.0;

	for(i = 0; i < 3; i++){
		distances_nucleus[0] += pow(positions[0][i],2);
		distances_nucleus[1] += pow(positions[1][i],2);
	}

	distances_nucleus[0] = sqrt(distances_nucleus[0]);
	distances_nucleus[1] = sqrt(distances_nucleus[1]);
}


// Function that caculated the auto-correlation function, the statistical inefficiency
void error_corr_func(double *A, int length, double *result){

	// Declaration and initiation of variables
	int i, k;
	double mean = 0;
	double mean2 = 0;
	double s = 0;
	double sigmaTot;
	int steps = 100;

	// Declaration of arrays
	double first_term[steps];
	double corr_func[steps];

	// Initiate the array first_term
	for(k = 0; k < steps; k++){
		first_term[k] = 0.0;
	}

	// Calculate all the expected values of A
	for(i = 0; i < length; i++){
		mean += A[i]/length;
		mean2 += ((A[i]*A[i])/length); 
	}

	// Calculate the first term
	for(i = 0; i < (length-steps); i++){
		for(k = 0; k < steps; k++){
			first_term[k] += (A[i]*A[i+k])/(length-steps);
		}
	}

	// Calculate the correlation function
	for(k = 0; k < steps; k++){
		corr_func[k] = ((first_term[k] - (mean*mean))/(mean2 - (mean*mean)));
	}

	// Calculate the statistical inefficiency
	i = 0;
	while(corr_func[i] >= exp(-2)){
		i++;
	}
	
	s = i;

	// Calculate the standard deviation
	sigmaTot = sqrt(s*(mean2 - mean*mean)/length);

	// Print the results in the terminal
	printf("Result: %.6f ± %.6f \n", mean, sigmaTot);
	printf("Statistical inefficiency (corr): %F \n", s);

	// Save the results in an array
	result[0] = mean;
	result[1] = sigmaTot;

}


// Calculate the statistical inefficiency from the block average
void error_block_average(double *A, int length){

	// Declaration and initiation of variables
	int i, j;
	int block_size;
	double mean, mean2, var_f;
	double mean_F, mean2_F, var_F;
	double s;

	// Create file to save data
	FILE *block;
	block = fopen("block_s.data","w");

	fprintf(block,"%f \n", 0);
	
	// Calculate statistical inefficiency for different block sizes
	for(block_size = 10; block_size < 500; block_size = block_size + 10){

		int nbr_blocks = length/block_size;
		double block_means[nbr_blocks];
		double block_means2[nbr_blocks];

		// Determine variance for the whole array
		mean = 0.0;
		mean2 = 0.0;	
		for(i = 0; i < length; i++){
			mean += A[i]/length;
			mean2 += A[i]*A[i]/length;
		}

		var_f = (mean2 - mean*mean);
		//printf("var: %.10f \n", var_f);
	
		// Determine average in each block
		for(i = 0; i < nbr_blocks; i++){
			
			block_means[i] = 0.0;
			block_means2[i] = 0.0;
			
			for(j = 0; j < block_size; j++){
				block_means[i] += A[(i*block_size+j)]/block_size;
			}
		}
	

		// Determine average in each block
		mean_F = 0.0;
		mean2_F = 0.0;
		for(i = 0; i < nbr_blocks; i++){
			mean_F += block_means[i]/nbr_blocks;
			mean2_F += block_means[i]*block_means[i]/nbr_blocks;
		}

		var_F = (mean2_F - mean_F * mean_F);
		s = block_size * var_F / var_f;
		
		fprintf(block,"%f \n", s);

	}

	// Print the statistical inefficiency in the terminal
	printf("Statistical inefficiency (block): %f \n", s);

	// Close file
	fclose(block);
}

// Function that returns the gradient of the ln(wavefunction) with respect to alpha
double get_grad_ln_wave(double distance, double alpha){
	
	double grad_ln_wave;

	grad_ln_wave = -distance * distance / pow(1 + alpha * distance, 2);

	return grad_ln_wave;

}

// Function that rescales the alpha-value
double rescale_alpha(double alpha, double energy_l[], double grad_ln_wave[], double distance, int iteration, int start_average, double beta){
	
	int i;
	double grad, gamma;
	double A = 1.0;
	double first_term = 0;
	double second_term = 0;
	double second_term_1 = 0;
	double second_term_2 = 0;
	double new_alpha;

	gamma = A * pow(iteration/start_average, -1*beta);

	for(i = iteration - start_average; i < iteration; i++){
		first_term += (energy_l[i] * grad_ln_wave[i]);
		second_term_1 += energy_l[i];
		second_term_2 += grad_ln_wave[i];
	}

	first_term = (double) first_term / start_average;
	second_term = (double) second_term_1 * second_term_2 / pow(start_average,2);

	grad = 2.0 * (first_term - second_term);

	new_alpha = alpha - (gamma * grad);

	return new_alpha;
}

