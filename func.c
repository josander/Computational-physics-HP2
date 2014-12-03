/*
func.c
Contains functions for homeproblem 2/b

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589

double getDistance(double positions[][3]){
	double distance = 0;
	int i;


	for(i = 0; i < 3; i++){
		distance += pow(positions[0][i] - positions[1][i], 2);
	}
		
	return (sqrt(distance));	


}

