#include <math.h>
#include <float.h>

#include "../headers/nrutil.h"

void calculate_Omega(int n_x, int n_y, int half_nx, int half_ny,
double *kx, double *ky, double Ceff[3][3][3][3], double *Omega11, 
double *Omega12, double*Omega21, double *Omega22){

int i5,i6;

double Omega_inv11, Omega_inv12, Omega_inv21, Omega_inv22;
double det_Omega_inv;
double n[3];
int J;

/* First step is to calculate the inverse of the acoustic tensor. Then
 * we invert it - We assume plane elasticity,  and hence matrix inversion
 * is straightforward */

for(i5=0; i5<n_x; ++i5){
	 n[1] = kx[i5];
for(i6=0; i6<n_y; ++i6){
	 n[2] = ky[i6];
		J = i6+n_y*i5;
		Omega_inv11 = Ceff[1][1][1][1]*n[1]*n[1]+
			Ceff[1][1][2][1]*n[1]*n[2]+
			Ceff[1][2][1][1]*n[2]*n[1]+
			Ceff[1][2][2][1]*n[2]*n[2];
		Omega_inv12 = Ceff[1][1][1][2]*n[1]*n[1]+
			Ceff[1][1][2][2]*n[1]*n[2]+
			Ceff[1][2][1][2]*n[2]*n[1]+
			Ceff[1][2][2][2]*n[2]*n[2];
		Omega_inv21 = Ceff[2][1][1][1]*n[1]*n[1]+
			Ceff[2][1][2][1]*n[1]*n[2]+
			Ceff[2][2][1][1]*n[2]*n[1]+
			Ceff[2][2][2][1]*n[2]*n[2];
		Omega_inv22 = Ceff[2][1][1][2]*n[1]*n[1]+
			Ceff[2][1][2][2]*n[1]*n[2]+
			Ceff[2][2][1][2]*n[2]*n[1]+
			Ceff[2][2][2][2]*n[2]*n[2];
			
	det_Omega_inv = Omega_inv11*Omega_inv22 - Omega_inv12*Omega_inv21;
	if(det_Omega_inv != 0.0){
		Omega11[J] = Omega_inv22/det_Omega_inv;  
		Omega22[J] = Omega_inv11/det_Omega_inv;  
		Omega12[J] = -Omega_inv12/det_Omega_inv;  
		Omega21[J] = -Omega_inv21/det_Omega_inv;  
	}
	else{
		Omega11[J] = Omega_inv22;
		Omega22[J] = Omega_inv11;
		Omega12[J] = -1.0*Omega_inv12;
		Omega21[J] = -1.0*Omega_inv21;
	}

}}

}
