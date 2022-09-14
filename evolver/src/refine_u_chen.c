/*********************************************************************
In an attempt to avoid unnecessary algebra, we have simplified the
summations to include only the non-zero terms while assuming that the
elastic constants are cubic and the eigenstrain is dilatational
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void refine_u_chen(int n_x, int n_y, int half_nx, int half_ny,
double delta_x, double delta_y, double **epsilon_T,
double delta_kx, double delta_ky, fftw_complex *comp,
double ave_comp, double Ceff[3][3][3][3], double DeltaC[3][3][3][3], 
double **sigma_zero, fftw_complex *u1_old, fftw_complex *u2_old, 
fftw_complex *u1_new, fftw_complex *u2_new, double **Del_sigma_T, 
double **sig_app, double **E, int MAXITS, double MAXERR,
int n_alpha, int n_beta, double *alpha, double *beta,
double *Omega11, double *Omega12, double *Omega21, double *Omega22,
fftw_plan planF, fftw_plan planB){

#ifndef DEBUG
FILE *fpw, *fp2;
#endif

FILE *fp1;

int i1,i2,i3,i4;

int ITS=0;
double ERR=1.0;
double tmp;

double *kx, *ky;

double zeta11, zeta12, zeta21, zeta22;

int nxny, J;
double inv_nxny;

fftw_complex gamma0, gamma1, gamma2, gamma3;

fftw_complex *u1_tmp, *u2_tmp;
fftw_complex *chi11, *chi12, *chi21, *chi22;
fftw_complex *Alpha, *Beta, *AlpBeta;
fftw_complex *eps_star11, *eps_star12, *eps_star22;

double a;

nxny = n_x*n_y;
inv_nxny = 1.0/nxny;
u1_tmp = fftw_malloc(nxny* sizeof(fftw_complex));
u2_tmp = fftw_malloc(nxny* sizeof(fftw_complex));
chi11 = fftw_malloc(nxny* sizeof(fftw_complex));
chi12 = fftw_malloc(nxny* sizeof(fftw_complex));
chi21 = fftw_malloc(nxny* sizeof(fftw_complex));
chi22 = fftw_malloc(nxny* sizeof(fftw_complex));
Alpha = fftw_malloc(nxny* sizeof(fftw_complex));
Beta = fftw_malloc(nxny* sizeof(fftw_complex));
AlpBeta = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star11 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star12 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star22 = fftw_malloc(nxny* sizeof(fftw_complex));

kx = (double *) malloc((size_t) n_x* sizeof(double));
ky = (double *) malloc((size_t) n_y* sizeof(double));

for(i1=0; i1 < n_x; ++i1){
if(i1 < half_nx) kx[i1] = i1*delta_kx;
else kx[i1] = (i1-n_x)*delta_kx;
}
for(i2=0; i2 < n_y; ++i2){
if(i2 < half_ny) ky[i2] = i2*delta_ky;
else ky[i2] = (i2-n_y)*delta_ky;
}


/* Calculate the derivative of displacement */

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
  chi11[J] = _Complex_I*kx[i1]*u1_old[J];
	chi12[J] = _Complex_I*ky[i2]*u1_old[J];
	chi21[J] = _Complex_I*kx[i1]*u2_old[J];
	chi22[J] = _Complex_I*ky[i2]*u2_old[J];
	eps_star12[J] = 0.5*(chi12[J]+chi21[J]);
}}

fftw_execute_dft(planB,chi11,chi11);
fftw_execute_dft(planB,chi12,chi12);
fftw_execute_dft(planB,chi21,chi21);
fftw_execute_dft(planB,chi22,chi22);
fftw_execute_dft(planB,eps_star12,eps_star12);

/* Multiply the derivative of displacement by alpha */

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){

	J = i2+n_y*i1;
	chi11[J] = chi11[J]*inv_nxny;
	chi12[J] = chi12[J]*inv_nxny;
	chi21[J] = chi21[J]*inv_nxny;
	chi22[J] = chi22[J]*inv_nxny;
	eps_star11[J] = chi11[J];
	eps_star12[J] = eps_star12[J]*inv_nxny;
	eps_star22[J] = chi22[J];

	a = alpha[(int) (n_alpha*creal(comp[J]))];


	chi11[J] = a*chi11[J];
	chi12[J] = a*chi12[J];
	chi21[J] = a*chi21[J];
	chi22[J] = a*chi22[J];
}}

fftw_execute_dft(planF,chi11,chi11);
fftw_execute_dft(planF,chi12,chi12);
fftw_execute_dft(planF,chi21,chi21);
fftw_execute_dft(planF,chi22,chi22);

/* Calculate the displacements, alpha, beta and alpha*beta */

fftw_execute_dft(planB,u1_old,u1_old);
fftw_execute_dft(planB,u2_old,u2_old);

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
	J = i2+n_y*i1;
	u1_old[J] = u1_old[J]*inv_nxny;
	u2_old[J] = u2_old[J]*inv_nxny;
	
	u1_tmp[J] = u1_old[J];
	u2_tmp[J] = u2_old[J];
	
	Alpha[J] = alpha[(int) (n_alpha*creal(comp[J]))];
	Beta[J] = beta[(int) (n_beta*creal(comp[J]))];
	AlpBeta[J] = Alpha[J]*Beta[J];
}}
fftw_execute_dft(planF,Alpha,Alpha);
fftw_execute_dft(planF,Beta,Beta);
fftw_execute_dft(planF,AlpBeta,AlpBeta);

zeta11 =  DeltaC[1][1][1][1]*E[1][1]
	+ DeltaC[1][1][2][2]*E[2][2];

zeta12 =  DeltaC[1][2][1][2]*E[1][2]
	+ DeltaC[1][2][2][1]*E[2][1];

zeta21 =  DeltaC[2][1][1][2]*E[1][2]
	+ DeltaC[2][1][2][1]*E[2][1];

zeta22 =  DeltaC[2][2][1][1]*E[1][1]
	+ DeltaC[2][2][2][2]*E[2][2];

/* Refine the solution */

while(ITS != MAXITS && ERR > MAXERR ){

	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
		
		J = i2+n_y*i1;

		gamma0 = DeltaC[1][1][1][1]*chi11[J]
							+ DeltaC[1][1][2][2]*chi22[J];
		gamma1 = DeltaC[1][2][1][2]*chi12[J]
							+ DeltaC[1][2][2][1]*chi21[J];
		gamma2 = DeltaC[2][1][1][2]*chi12[J]
							+ DeltaC[2][1][2][1]*chi21[J];
		gamma3 = DeltaC[2][2][1][1]*chi11[J]
							+ DeltaC[2][2][2][2]*chi22[J];

		u1_new[J] =  - _Complex_I*Omega11[J]*kx[i1]
		*(sigma_zero[1][1]*Beta[J]
		-zeta11*Alpha[J]
		+ Del_sigma_T[1][1]*AlpBeta[J] - gamma0)
		- _Complex_I*Omega21[J]*kx[i1]
		*(sigma_zero[2][1]*Beta[J]
		-zeta21*Alpha[J]
		+ Del_sigma_T[2][1]*AlpBeta[J] - gamma2)
		- _Complex_I*Omega11[J]*ky[i2]
		*(sigma_zero[1][2]*Beta[J]
		-zeta12*Alpha[J]
		+ Del_sigma_T[1][2]*AlpBeta[J] - gamma1)
		- _Complex_I*Omega21[J]*ky[i2]
		*(sigma_zero[2][2]*Beta[J]
		-zeta22*Alpha[J]
		+ Del_sigma_T[2][2]*AlpBeta[J] - gamma3);
		u2_new[J] =  - _Complex_I*Omega12[J]*kx[i1]
		*(sigma_zero[1][1]*Beta[J]
		-zeta11*Alpha[J]
		+ Del_sigma_T[1][1]*AlpBeta[J] - gamma0)
		- _Complex_I*Omega22[J]*kx[i1]
		*(sigma_zero[2][1]*Beta[J]
		-zeta21*Alpha[J]
		+ Del_sigma_T[2][1]*AlpBeta[J] - gamma2)
		- _Complex_I*Omega12[J]*ky[i2]
		*(sigma_zero[1][2]*Beta[J]
		-zeta12*Alpha[J]
		+ Del_sigma_T[1][2]*AlpBeta[J] - gamma1)
		- _Complex_I*Omega22[J]*ky[i2]
		*(sigma_zero[2][2]*Beta[J]
		-zeta22*Alpha[J]
		+ Del_sigma_T[2][2]*AlpBeta[J] - gamma3);
		
	chi11[J] = _Complex_I*kx[i1]*u1_new[J];
	chi12[J] = _Complex_I*ky[i2]*u1_new[J];
	chi21[J] = _Complex_I*kx[i1]*u2_new[J];
	chi22[J] = _Complex_I*ky[i2]*u2_new[J];
	eps_star12[J] = 0.5*(chi12[J]+chi21[J]);
	
	}}

	fftw_execute_dft(planB,chi11,chi11);
	fftw_execute_dft(planB,chi12,chi12);
	fftw_execute_dft(planB,chi21,chi21);
	fftw_execute_dft(planB,chi22,chi22);
	fftw_execute_dft(planB,eps_star12,eps_star12);
	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
	J = i2+n_y*i1;
	chi11[J] = chi11[J]*inv_nxny;  
	chi12[J] = chi12[J]*inv_nxny;
	chi21[J] = chi21[J]*inv_nxny;
	chi22[J] = chi22[J]*inv_nxny;
	eps_star12[J] = eps_star12[J]*inv_nxny;
	eps_star11[J] = chi11[J];
	eps_star22[J] = chi22[J];
	a = alpha[(int) (n_alpha*creal(comp[J]))];
	chi11[J] = a*chi11[J];  
	chi12[J] = a*chi12[J];
	chi21[J] = a*chi21[J];
	chi22[J] = a*chi22[J];
	}}
	fftw_execute_dft(planF,chi11,chi11);
	fftw_execute_dft(planF,chi12,chi12);
	fftw_execute_dft(planF,chi21,chi21);
	fftw_execute_dft(planF,chi22,chi22);

	for(i3=0; i3<n_x; ++i3){
	for(i4=0; i4<n_y; ++i4){
		J = i4+n_y*i3;
		u1_tmp[J] = u1_old[J];
		u2_tmp[J] = u2_old[J];
		u1_old[J] = u1_new[J];
		u2_old[J] = u2_new[J];
	}}
	
	fftw_execute_dft(planB,u1_old,u1_old);
	fftw_execute_dft(planB,u2_old,u2_old);
	
	ERR = 0.0;
	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
		J = i2+n_y*i1;	
		u1_old[J] = u1_old[J]*inv_nxny;  
		u2_old[J] = u2_old[J]*inv_nxny; 
		tmp  = creal((u1_old[J]-u1_tmp[J])
								*(u1_old[J]-u1_tmp[J])
								+(u2_old[J]-u2_tmp[J])
								*(u2_old[J]-u2_tmp[J]));
		ERR = ERR + tmp;
	}}

/* Calculate the error */
	
	ERR = sqrt(ERR);
	ERR = ERR*inv_nxny;
	if(ERR > 1.0){
	printf("The ERR is %le\n",ERR);
	printf("It cannot be more than unity. Hence exiting\n");
	exit(0);
	}
	
	ITS = ITS + 1;
}

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
	J = i2+n_y*i1;
	u1_old[J] = u1_new[J];
	u2_old[J] = u2_new[J];
}}

#ifndef DEBUG
fpw = fopen("../output/itnos_err","a");
fprintf(fpw,"%d %le\n",ITS-1,ERR);
fclose(fpw);
if(ITS == MAXITS){
	printf("No convergence after %d iterations\n",MAXITS);
	printf("Exiting.\n"); 
	printf("MAXITS:\n%d Error:%le\n",MAXITS,ERR);
	fp1 = fopen("../output/unconverged/README","w");
	fprintf(fp1,"MAXITS:\n%d\n",MAXITS);
  fprintf(fp1,"n_x:\n%d\nn_y:\n%d\n",n_x,n_y);
  fclose(fp1);
  fp1 = fopen("../output/unconverged/u1","w");
  fp2 = fopen("../output/unconverged/u2","w");
  for(i1=0; i1<n_x; ++i1){
  for(i2=0; i2<n_y; ++i2){
    J = i2+n_y*i1;
    fprintf(fp1,"%le %le\n",creal(u1_new[J]),cimag(u1_new[J]));
	  fprintf(fp2,"%le %le\n",creal(u2_new[J]),cimag(u2_new[J]));
  }}
  fclose(fp1);
  fclose(fp2);
	exit(0);
}
#endif

#ifdef DEBUG
fp1 = fopen("output/itnos_err.dat","w");
fprintf(fp1,"%d %le\n",ITS-1,ERR);
fclose(fp1);
#endif

free(kx);
free(ky);
fftw_free(u1_tmp);
fftw_free(u2_tmp);
fftw_free(chi11);
fftw_free(chi12);
fftw_free(chi21);
fftw_free(chi22);
fftw_free(Alpha);
fftw_free(Beta);
fftw_free(AlpBeta);
fftw_free(eps_star11);
fftw_free(eps_star12);
fftw_free(eps_star22);

}
