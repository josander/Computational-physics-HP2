/*
func.h
*/

#ifndef _func_h
#define _func_h

extern double getDistance(double [][3]);
extern double get_local_e(double [][3], double);
extern double get_wavefunction(double [][3], double, double);
extern void get_distances_nucleus(double [][3], double []);
extern void error_corr_func(double *, int);
extern void error_block_average(double *, int);
extern double rescale_alpha(double, double [], double [], double, int);
extern double get_grad_ln_wave(double, double);

#endif
