#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void evolve_hom_unconstrained(int n_x, int n_y, double delta_x, 
double delta_y, double kappa, double A, double delta_t, 
int time_steps, double t0, fftw_complex *comp, double Ceff[3][3][3][3], 
double **sigma_T, double **epsilon_T, double ave_comp, double **sig_app, 
int STEPS, int n_beta, double *beta, double *beta_prime, 
int plan_indicator, int SNA, int NOISEADDER, double noise_str){

FILE *fpd;
FILE *fptmp;

int INDEX=0;
int i1,i2,i3,i4,i5,i6;
int J;

int half_nx, half_ny;
double *kx, *ky;
double delta_kx, delta_ky;
double k2, k4;
double inv_denom;

int nxny;
double inv_nxny;

char NAME[50];

double *c;
double **E;
double **strain;
double *Omega11;
double *Omega12;
double *Omega21;
double *Omega22;
double vol_fraction;
double S[3][3][3][3];

double b;
double bp;
double realc;
size_t tmp;

unsigned FLAG;

fftw_complex *g;
fftw_complex *u1,*u2;
fftw_complex *eps_star11, *eps_star12, *eps_star22;
fftw_complex *eta;
fftw_complex *mu_el;

fftw_plan planF, planB;

nxny = n_x*n_y;
inv_nxny = 1.0/nxny;

g = fftw_malloc(nxny* sizeof(fftw_complex));
u1 = fftw_malloc(nxny* sizeof(fftw_complex));
u2 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star11 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star12 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star22 = fftw_malloc(nxny* sizeof(fftw_complex));
eta = fftw_malloc(2*2* sizeof(fftw_complex));
mu_el = fftw_malloc(nxny* sizeof(fftw_complex));

FLAG = FFTW_EXHAUSTIVE;
if(plan_indicator == 0){
FLAG = FFTW_ESTIMATE;
}
else if(plan_indicator == 1){
FLAG = FFTW_MEASURE;
}
else if(plan_indicator == 2){
FLAG = FFTW_PATIENT;
}
else if(plan_indicator == 3){
FLAG = FFTW_EXHAUSTIVE;
}

planF =
fftw_plan_dft_2d(n_x,n_y,mu_el,mu_el,FFTW_FORWARD,FLAG);
planB =
fftw_plan_dft_2d(n_x,n_y,eps_star11,eps_star11,FFTW_BACKWARD,FLAG);

/* Declare and initialise the acoustic tensor */

Omega11 = dvector(0,nxny);
Omega12 = dvector(0,nxny);
Omega21 = dvector(0,nxny);
Omega22 = dvector(0,nxny);

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
J = i2+n_y*i1;
Omega11[J] = 0.0;
Omega12[J] = 0.0;
Omega21[J] = 0.0;
Omega22[J] = 0.0;
}}

/* Write the initial configuration */

c = (double *) malloc((size_t) nxny* sizeof(double));

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
J = i2+n_y*i1;
c[J] = __real__ comp[J];
}}

sprintf(NAME,"../output/data/time%d.dat", (int) (INDEX + t0));
fpd = fopen(NAME,"w");
tmp = fwrite(&c[0],sizeof(double),(size_t) nxny,fpd);
fclose(fpd);

/* Calculate homogeneous strain */

strain = dmatrix(1,2,1,2);
E = dmatrix(1,2,1,2);

calculate_S(Ceff,S);

vol_fraction = 0.0;
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2+n_y*i1;
realc = __real__ comp[J];
b = beta[(int) (n_beta*realc)];
vol_fraction = vol_fraction + b;
}}
vol_fraction = vol_fraction*inv_nxny;

for(i1=1; i1<3; ++i1){
for(i2=1; i2<3; ++i2){
E[i1][i2] = epsilon_T[i1][i2]*vol_fraction;
}}

for(i1=1; i1<3; ++i1){
for(i2=1; i2<3; ++i2){
for(i3=1; i3<3; ++i3){
for(i4=1; i4<3; ++i4){
E[i1][i2] = E[i1][i2] + S[i1][i2][i3][i4]*sig_app[i3][i4];
}}}}


/* Evolve */

half_nx = (int) n_x/2;
half_ny = (int) n_y/2;

delta_kx = (2.0*M_PI)/(n_x*delta_x);
delta_ky = (2.0*M_PI)/(n_y*delta_y);

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

/* Calculate g and its Fourier transform */

fptmp = fopen("../output/noise_info","a");

for(INDEX=0; INDEX < time_steps+1; ++INDEX){

printf("%d\n",INDEX);

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2+n_y*i1;
realc = __real__ comp[J];
g[J] = 2.0*A*realc*(1.0-realc)*(1.0-2.0*realc) + 0.0*_Complex_I;
}}
fftw_execute_dft(planF,g,g);

if(INDEX == 0){

/* The acoustic tensor Omega and the (homogeneous) base solution */

calculate_Omega(n_x,n_y,half_nx,half_ny,kx,ky,Ceff,
Omega11,Omega12,Omega21,Omega22);

}

calculate_uzero(n_x,n_y,half_nx,half_ny,delta_kx,delta_ky,ave_comp,
comp,Ceff,sigma_T,u1,u2,n_beta,beta,Omega11,Omega12,Omega21,
Omega22,planF);

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	eps_star11[J] = _Complex_I*kx[i1]*u1[J];
	eps_star22[J] = _Complex_I*ky[i2]*u2[J];
	eps_star12[J] = 0.5*_Complex_I*(kx[i1]*u2[J]+ky[i2]*u1[J]);
}}

/* Get the strain back to the real space */

fftw_execute(planB);
fftw_execute_dft(planB,eps_star12,eps_star12);
fftw_execute_dft(planB,eps_star22,eps_star22);
for(i3=0; i3<n_x; ++i3){
for(i4=0; i4<n_y; ++i4){
J = i4+n_y*i3;
eps_star11[J] = eps_star11[J]*inv_nxny;
eps_star22[J] = eps_star22[J]*inv_nxny;
eps_star12[J] = eps_star12[J]*inv_nxny;
}}

/* Calculate mu_el and take it to the Fourier space */

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	realc = __real__ comp[J];
	mu_el[J] = 0.0;
	
	b = beta[(int) (n_beta*realc)];
	bp = beta_prime[(int) (n_beta*realc)];
	
	strain[1][1] = 
	E[1][1] + __real__ eps_star11[J] - epsilon_T[1][1]*b;
	strain[1][2] = 
	E[1][2] + __real__ eps_star12[J] - epsilon_T[1][2]*b;
	strain[2][1] = 
	E[2][1] + __real__ eps_star12[J] - epsilon_T[2][1]*b;
	strain[2][2] = 
	E[2][2] + __real__ eps_star22[J] - epsilon_T[2][2]*b;

	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
	for(i5=1; i5 < 3; ++i5){
	for(i6=1; i6 < 3; ++i6){
	mu_el[J] = mu_el[J]
		-bp*epsilon_T[i3][i4]*strain[i5][i6]*Ceff[i3][i4][i5][i6];
	}}}}
}}

fftw_execute(planF);
fftw_execute_dft(planF,comp,comp);

/* Evolve the composition profile */


for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){

J = i2+n_y*i1;

k2 = kx[i1]*kx[i1] + ky[i2]*ky[i2];
k4 = k2*k2;

inv_denom = 1.0 + 2.0*kappa*k4*delta_t;
inv_denom = 1.0/inv_denom;

comp[J] = inv_denom*( comp[J] - k2*delta_t*(g[J]+mu_el[J]) );

}}

/* Get the composition back to real space */

fftw_execute_dft(planB,comp,comp);

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2+n_y*i1;
comp[J] = comp[J]*inv_nxny + 0.0*_Complex_I;
}}

/* Once in a while check that the composition is within bounds,
 * and write the output **/

if(INDEX != 0 && INDEX%STEPS ==0){
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
		J = i2+n_y*i1;
		c[J] = __real__ comp[J];
		if(c[J] <= -0.1 || c[J] >= 1.1){
		printf("The composition goes out of bounds\n");
		printf("Exiting\n");
		exit(0);
		}
	}}
	
	sprintf(NAME,"../output/data/time%d.dat",
		(int) (INDEX + t0));
	fpd = fopen(NAME,"w");
	tmp = fwrite(&c[0],sizeof(double),(size_t) nxny,fpd);
	fclose(fpd);

}

if(INDEX != 0 && INDEX%NOISEADDER == 0 && SNA == 1){
	fprintf(fptmp,"Step %d: Noise added\n",INDEX);	
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	c[J] = __real__ comp[J];
	}}
	add_noise(n_x,n_y,c,noise_str);
}

}
fclose(fptmp);

/* Free the variables and destroy the plans */

fftw_free(g);
fftw_free(eps_star11);
fftw_free(eps_star12);
fftw_free(eps_star22);
fftw_free(u1);
fftw_free(u2);
fftw_free(mu_el);
fftw_free(eta);

fftw_destroy_plan(planF);
fftw_destroy_plan(planB);

free(c);
free(kx);
free(ky);
free_dmatrix(E,1,2,1,2);
free_dmatrix(strain,1,2,1,2);
free_dvector(Omega11,0,nxny);
free_dvector(Omega12,0,nxny);
free_dvector(Omega21,0,nxny);
free_dvector(Omega22,0,nxny);
}
