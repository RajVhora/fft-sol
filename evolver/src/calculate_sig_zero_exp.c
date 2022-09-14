/*************************************************************
We previously assumed that the eigenstrain is dilatational and the elastic
constants are cubic, and hence shortened the lengthy summations, 
including only the non-zero terms. However, now the sig_zero
expectation values are calculated without making any such assumption
and hence skipping terms.
M P Gururajan and N Davessar, 1 April, 20
***************************************************************/

#include<complex.h>
#include<fftw3.h>

#include "../headers/functions.h"

void calculate_sig_zero_exp(int n_x, int n_y, double delta_x, 
double delta_y, fftw_complex *comp, double ave_comp, 
double Ceff[3][3][3][3], double DeltaC[3][3][3][3], double **epsilon_T,
double **sig_zero_exp, int n_alpha, int n_beta, double *alpha, 
double *beta){

int i5,i6;
double c;
double a, b, d, e;

/* Calculate the 11, 22 and 12 component of 'expectation' value of eigenstrain */

sig_zero_exp[1][1] = 0.0;
sig_zero_exp[2][2] = 0.0;
sig_zero_exp[1][2] = 0.0;
for(i5=0; i5 < n_x; ++i5){
for(i6=0; i6 < n_y; ++i6){
c = creal(comp[i6+n_y*i5]);
a = alpha[(int) (n_alpha*c)];
b = beta[(int) (n_beta*c)] * epsilon_T[1][1];
d = beta[(int) (n_beta*c)] * epsilon_T[2][2];
e = beta[(int) (n_beta*c)] * epsilon_T[1][2];
sig_zero_exp[1][1] = sig_zero_exp[1][1] 
+ (Ceff[1][1][1][1]+DeltaC[1][1][1][1]*a)*b 
+ (Ceff[1][1][2][2]+DeltaC[1][1][2][2]*a)*d 
+ (Ceff[1][1][1][2]+DeltaC[1][1][1][2]*a)*e
+ (Ceff[1][1][2][1]+DeltaC[1][1][2][1]*a)*e;

sig_zero_exp[2][2] = sig_zero_exp[2][2] 
+ (Ceff[2][2][1][1]+DeltaC[2][2][1][1]*a)*b 
+ (Ceff[2][2][2][2]+DeltaC[2][2][2][2]*a)*d 
+ (Ceff[2][2][1][2]+DeltaC[2][2][1][2]*a)*e
+ (Ceff[2][2][2][1]+DeltaC[2][2][2][1]*a)*e;

sig_zero_exp[1][2] = sig_zero_exp[1][2] 
+ (Ceff[1][2][1][1]+DeltaC[1][2][1][1]*a)*b 
+ (Ceff[1][2][2][2]+DeltaC[1][2][2][2]*a)*d 
+ (Ceff[1][2][1][2]+DeltaC[1][2][1][2]*a)*e
+ (Ceff[1][2][2][1]+DeltaC[1][2][2][1]*a)*e;
}}

sig_zero_exp[1][1] = sig_zero_exp[1][1]/(n_x*n_y);
sig_zero_exp[2][2] = sig_zero_exp[2][2]/(n_x*n_y);
sig_zero_exp[1][2] = sig_zero_exp[1][2]/(n_x*n_y);
sig_zero_exp[2][1] = sig_zero_exp[1][2];

}
