#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void evolve_inhom_unconst_varmob(int n_x, int n_y, double delta_x, 
double delta_y, double kappa, double A, double M_b, double M_n, 
double M_t, double delta_t, int time_steps, 
double t0, fftw_complex *comp, double Ceff[3][3][3][3], 
double DeltaC[3][3][3][3], double **sigma_T, double **epsilon_T, 
double ave_comp, double **sig_app, int STEPS, int MAXITR, double MAXERR, 
int n_alpha, int n_beta, double *alpha, double *beta, double *alpha_prime, 
double *beta_prime, int plan_indicator, int SNA, int NOISEADDER, 
double noise_str){

FILE *fpd;
FILE *fptmp;

int INDEX=0;
int i1,i2,i3,i4,i5,i6;
int J;

int half_nx, half_ny;
double *kx, *ky;
double delta_kx, delta_ky;
double k2, k4;
double denom,inv_denom;
double M_ave;
double M[3][3], n[3];
double phi_c, mod_n;

int nxny;
double inv_nxny;

char NAME[50];

double epsstr11, epsstr22, epsstr12;

double *c;
double **Del_sigma_T;
double **E;
double **strain;
double *Omega11;
double *Omega12;
double *Omega21;
double *Omega22;

double a, b;
double ap, bp;
double realc;
size_t tmp;

unsigned FLAG;

fftw_complex *g;
fftw_complex *u1_old,*u2_old;
fftw_complex *u1_new,*u2_new;
fftw_complex *eps_star11, *eps_star12, *eps_star22;
fftw_complex *eta;
fftw_complex *mu_el;
fftw_complex *delc1,*delc2,*h1,*h2,*Mh1,*Mh2;

fftw_plan planF, planB;

nxny = n_x*n_y;
inv_nxny = 1.0/nxny;

M_ave = 0.5*(M_b + M_t);

g = fftw_malloc(nxny* sizeof(fftw_complex));
u1_old = fftw_malloc(nxny* sizeof(fftw_complex));
u2_old = fftw_malloc(nxny* sizeof(fftw_complex));
u1_new = fftw_malloc(nxny* sizeof(fftw_complex));
u2_new = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star11 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star12 = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star22 = fftw_malloc(nxny* sizeof(fftw_complex));
eta = fftw_malloc(2*2* sizeof(fftw_complex));
mu_el = fftw_malloc(nxny* sizeof(fftw_complex));
delc1 = fftw_malloc(nxny* sizeof(fftw_complex));
delc2 = fftw_malloc(nxny* sizeof(fftw_complex));
h1 = fftw_malloc(nxny* sizeof(fftw_complex));
h2 = fftw_malloc(nxny* sizeof(fftw_complex));
Mh1 = fftw_malloc(nxny* sizeof(fftw_complex));
Mh2 = fftw_malloc(nxny* sizeof(fftw_complex));

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

half_nx = (int) n_x/2;
half_ny = (int) n_y/2;

kx = (double *) malloc((size_t) n_x* sizeof(double));
ky = (double *) malloc((size_t) n_y* sizeof(double));

delta_kx = (2.0*M_PI)/(n_x*delta_x);
delta_ky = (2.0*M_PI)/(n_y*delta_y);

Omega11 = dvector(0,nxny);
Omega12 = dvector(0,nxny);
Omega21 = dvector(0,nxny);
Omega22 = dvector(0,nxny);

for(i1=0; i1 < n_x; ++i1){
if(i1 < half_nx) kx[i1] = i1*delta_kx;
else kx[i1] = (i1-n_x)*delta_kx;
}
for(i2=0; i2 < n_y; ++i2){
if(i2 < half_ny) ky[i2] = i2*delta_ky;
else ky[i2] = (i2-n_y)*delta_ky;
}

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
c[J] = creal(comp[J]);
}}

sprintf(NAME,"../output/data/time%d.dat",
(int) (INDEX + t0));
fpd = fopen(NAME,"w");
tmp = fwrite(&c[0],sizeof(double),(size_t) nxny,fpd);
fclose(fpd);

/** Calculate the Del_sigma_T tensor **/

Del_sigma_T = dmatrix(1,2,1,2);
calculate_Del_sigma_T(DeltaC,epsilon_T,Del_sigma_T);
E = dmatrix(1,2,1,2);
strain = dmatrix(1,2,1,2);

/* Evolve */


/* Calculate g and its Fourier transform */

fptmp = fopen("../output/noise_info","a");

for(INDEX=0; INDEX < time_steps+1; ++INDEX){

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2+n_y*i1;
realc = creal(comp[J]);
g[J] = 2.0*A*realc*(1.0-realc)*(1.0-2.0*realc) + 0.0*_Complex_I;
}}
fftw_execute_dft(planF,g,g);

if(INDEX == 0){

/* The acoustic tensor Omega and the (homogeneous) base solution */

calculate_Omega(n_x,n_y,half_nx,half_ny,kx,ky,Ceff,
Omega11,Omega12,Omega21,Omega22);

calculate_uzero(n_x,n_y,half_nx,half_ny,delta_kx,delta_ky,ave_comp,
comp,Ceff,sigma_T,u1_old,u2_old,n_beta,beta,Omega11,Omega12,Omega21,
Omega22,planF);

}


/* Refine the displacements */

refine_u(n_x,n_y,half_nx,half_ny,delta_x,delta_y,epsilon_T,
delta_kx,delta_ky,comp,ave_comp,Ceff,DeltaC,sigma_T,
u1_old,u2_old,u1_new,u2_new,Del_sigma_T,sig_app,E,MAXITR,MAXERR,
n_alpha,n_beta,alpha,beta,Omega11,Omega12,Omega21,Omega22,planF,planB);


for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	eps_star11[J] = _Complex_I*kx[i1]*u1_new[J];
	eps_star22[J] = _Complex_I*ky[i2]*u2_new[J];
	eps_star12[J] = 0.5*_Complex_I*(kx[i1]*u2_new[J]
																		+ky[i2]*u1_new[J]);
}}

/* Get the strain back to the real space */

fftw_execute(planB);
fftw_execute_dft(planB,eps_star12,eps_star12);
fftw_execute_dft(planB,eps_star22,eps_star22);

/* Calculate mu_el and take it to the Fourier space */

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){

	J = i2+n_y*i1;
	
	epsstr11 = creal(eps_star11[J]*inv_nxny);
	epsstr22 = creal(eps_star22[J]*inv_nxny);
	epsstr12 = creal(eps_star12[J]*inv_nxny);
	
	realc = creal(comp[J]);
	c[J] = realc;
	mu_el[J] = 0.0;
	
	a = alpha[(int) (n_alpha*realc)];
	b = beta[(int) (n_beta*realc)];
	
	ap = alpha_prime[(int) (n_alpha*realc)];
	bp = beta_prime[(int) (n_beta*realc)];
	
	strain[1][1] = 
	E[1][1] + epsstr11 - epsilon_T[1][1]*b;
	strain[1][2] = 
	E[1][2] + epsstr12 - epsilon_T[1][2]*b;
	strain[2][1] = 
	E[2][1] + epsstr12 - epsilon_T[2][1]*b;
	strain[2][2] = 
	E[2][2] + epsstr22 - epsilon_T[2][2]*b;

	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
	for(i5=1; i5 < 3; ++i5){
	for(i6=1; i6 < 3; ++i6){
	mu_el[J] = mu_el[J]
		+0.5*ap*DeltaC[i3][i4][i5][i6]*strain[i5][i6]*strain[i3][i4]
		-bp*epsilon_T[i3][i4]*strain[i5][i6]
			*(Ceff[i3][i4][i5][i6]+DeltaC[i3][i4][i5][i6]*a);
	}}}}
}}

fftw_execute(planF);
fftw_execute_dft(planF,comp,comp);

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	delc1[J] = _Complex_I*kx[i1]*comp[J];
	delc2[J] = _Complex_I*ky[i2]*comp[J];
}}
fftw_execute_dft(planB,delc1,delc1);
fftw_execute_dft(planB,delc2,delc2);
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	delc1[J] = creal(delc1[J]*inv_nxny);
	delc2[J] = creal(delc2[J]*inv_nxny);
}}

/* Evolve the composition profile */

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
J = i2+n_y*i1;
k2 = kx[i1]*kx[i1] + ky[i2]*ky[i2];
h1[J] = -_Complex_I*kx[i1]*(g[J]+2.0*kappa*k2*comp[J]+mu_el[J]);
h2[J] = -_Complex_I*ky[i2]*(g[J]+2.0*kappa*k2*comp[J]+mu_el[J]);
}}
fftw_execute_dft(planB,h1,h1);
fftw_execute_dft(planB,h2,h2);
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	h1[J] = creal(h1[J]*inv_nxny);
	h2[J] = creal(h2[J]*inv_nxny);
}}

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
J = i2+n_y*i1;
phi_c = c[J]*c[J]*(1.0-c[J])*(1.0-c[J]);
mod_n = delc1[J]*delc1[J] + delc2[J]*delc2[J];
mod_n = sqrt(mod_n);
if(mod_n > 1.0e-6){
n[1] = delc1[J]/mod_n;
n[2] = delc2[J]/mod_n;
}
else{
n[1] = 0.0;
n[2] = 0.0;
}
M[1][1] = M_b + phi_c*(M_t*n[2]*n[2] + M_n*n[1]*n[1]);
M[2][2] = M_b + phi_c*(M_t*n[1]*n[1] + M_n*n[2]*n[2]);
M[1][2] = phi_c*n[1]*n[2]*(M_n-M_t);
M[2][1] = M[1][2];
Mh1[J] = M[1][1]*h1[J] + M[1][2]*h2[J];
Mh2[J] = M[2][1]*h1[J] + M[2][2]*h2[J];
}}
fftw_execute_dft(planF,Mh1,Mh1);
fftw_execute_dft(planF,Mh2,Mh2);

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
J = i2+n_y*i1;
k2 = kx[i1]*kx[i1] + ky[i2]*ky[i2];
k4 = k2*k2;
denom = 1.0 + M_ave*kappa*k4*delta_t;
inv_denom = 1.0/denom;

comp[J] = inv_denom*
	(comp[J]*denom-_Complex_I*delta_t*(Mh1[J]*kx[i1]+Mh2[J]*ky[i2]));

}}

/* Get the composition back to real space */

fftw_execute_dft(planB,comp,comp);
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	comp[J] = comp[J]*inv_nxny;
}}

/* Once in a while check that the composition is within bounds,
 * and write the output **/

if(INDEX != 0 && INDEX%STEPS ==0){
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
		J = i2+n_y*i1;
		c[J] = creal(comp[J]);
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
	c[J] = creal(comp[J]);
	}}
	add_noise(n_x,n_y,c,noise_str);
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	comp[J] = c[J] + _Complex_I*0.0;
	}}
}

}
fclose(fptmp);

/* Free the variables and destroy the plans */

fftw_free(g);
fftw_free(eps_star11);
fftw_free(eps_star12);
fftw_free(eps_star22);
fftw_free(u1_old);
fftw_free(u2_old);
fftw_free(u1_new);
fftw_free(u2_new);
fftw_free(mu_el);
fftw_free(eta);
fftw_free(delc1);
fftw_free(delc2);
fftw_free(h1);
fftw_free(h2);
fftw_free(Mh1);
fftw_free(Mh2);

fftw_destroy_plan(planF);
fftw_destroy_plan(planB);

free(c);
free(kx);
free(ky);
free_dmatrix(Del_sigma_T,1,2,1,2);
free_dmatrix(E,1,2,1,2);
free_dmatrix(strain,1,2,1,2);
free_dvector(Omega11,0,nxny);
free_dvector(Omega12,0,nxny);
free_dvector(Omega21,0,nxny);
free_dvector(Omega22,0,nxny);
}
