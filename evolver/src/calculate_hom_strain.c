/*********************************************************************
We assume that the eigenstrain is dilatational and the elastic constants 
are cubic, and hence shortened the lengthy summations, including only 
the non-zero terms 
**********************************************************************/
/*********************************************************************
UPDATE:
The summation equation has been modified to include results for 
non-dilatational eigenstrains as well. Refer to values of A11 and A22
in line 56/75 and 58/59 respectively.
Modified by Nitin Davessar
01/04/2013
**********************************************************************/

#include<stdio.h>
#include<stdlib.h>

#include<complex.h>
#include<fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void calculate_hom_strain(int n_x, int n_y, 
double delta_x, double delta_y, double Ceff[3][3][3][3], 
double DeltaC[3][3][3][3], double S[3][3][3][3], 
double **sig_zero_exp, fftw_complex *comp, double ave_comp,
fftw_complex *eps_star11, fftw_complex *eps_star12,
fftw_complex *eps_star22, double **sig_app, double **E, 
int n_alpha, double *alpha){

int i1,i2;
int J;

double A11=0.0;
double A22=0.0;
double A12=0.0;
double A21=0.0;
double epsstr11,epsstr22,epsstr12;

double a;
double inv_nxny;

inv_nxny = 1.0/(n_x*n_y);

/* Calculate the expectation value for the periodic strain */

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2+n_y*i1;
epsstr11 = creal(eps_star11[J]);
epsstr22 = creal(eps_star22[J]);
epsstr12 = creal(eps_star12[J]);
a = alpha[(int) (n_alpha*creal(comp[J]))];

A11 = A11 + (Ceff[1][1][1][1]+DeltaC[1][1][1][1]*a)*epsstr11
+ (Ceff[1][1][2][2]+DeltaC[1][1][2][2]*a)*epsstr22+(Ceff[1][1][1][2]+DeltaC[1][1][1][2]*a)*epsstr12+(Ceff[1][1][2][1]+DeltaC[1][1][2][1]*a)*epsstr12;
A22 = A22 + (Ceff[2][2][1][1]+DeltaC[2][2][1][1]*a)*epsstr11
+ (Ceff[2][2][2][2]+DeltaC[2][2][2][2]*a)*epsstr22+(Ceff[2][2][1][2]+DeltaC[2][2][1][2]*a)*epsstr12+(Ceff[2][2][2][1]+DeltaC[2][2][2][1]*a)*epsstr12;
A12 = A12 + (Ceff[1][2][1][2]+DeltaC[1][2][1][2]*a
+ Ceff[1][2][2][1]+DeltaC[1][2][2][1]*a)*epsstr12;

}}

/* Calculate the 'effective' stress */

A11 = sig_zero_exp[1][1] + sig_app[1][1] - A11*inv_nxny;
A22 = sig_zero_exp[2][2] + sig_app[2][2] - A22*inv_nxny;
A12 = sig_app[1][2] - A12*inv_nxny;
A21 = A12;

/* Multiply the 'effective' stress by the 'effective' compliance to
 * obtain the homogeneous strain */

E[1][1] = S[1][1][1][1]*A11 + S[1][1][2][2]*A22;
E[1][2] = (S[1][2][1][2] + S[1][2][2][1])*A12;
E[2][1] = E[1][2];
E[2][2] = S[2][2][1][1]*A11 + S[2][2][2][2]*A22;

}
