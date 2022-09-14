#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void evolve_inhom_unconstrained_ceta(int n_x, int n_y, double delta_x, 
double delta_y, double kappac, double kappae, double A, double B, double P,
double M, double L,
double delta_t, int time_steps, 
double t0, fftw_complex *comp, fftw_complex *eta, double Ceff[3][3][3][3], 
double DeltaC[3][3][3][3], double **sigma_Tc, double **epsilon_Tc,
double **sigma_Te, double **epsilon_Te, 
double ave_comp, double **sig_app, int STEPS, int MAXITR, double MAXERR, 
int n_alpha, int n_beta, double *alpha, double *beta, double *alpha_prime, 
double *beta_prime, int plan_indicator, int SNA, int NOISEADDER, 
double noise_str){


FILE *fpd;
FILE *fpmu;
FILE *fptmp;
//
FILE *fit;
//
int INDEX=0;
int i1,i2,i3,i4,i5,i6;
int J;

int half_nx, half_ny;
double *kx, *ky;
double delta_kx, delta_ky;
double k2, k4;
double inv_denom;
//
double rr;
double s11, s12, s21, s22;
int i11, i22;
//
int nxny;
double inv_nxny;

char NAME[50];

double epsstr11c, epsstr22c, epsstr12c;
double epsstr11e, epsstr22e, epsstr12e;

double *c;
double *e;
double **Del_sigma_Tc;
double **Del_sigma_Te;
double **Ec;
double **strainc;
double *Omega11;
double *Omega12;
double *Omega21;
double *Omega22;
double **Ee;
double **straine;

double a, b;
double ap, bp;
double realc;
double reale;
size_t tmp;

double W, Wp;

unsigned FLAG;

fftw_complex *g;
fftw_complex *u1_oldc,*u2_oldc;
fftw_complex *u1_newc,*u2_newc;
fftw_complex *eps_star11c, *eps_star12c, *eps_star22c;
fftw_complex *mu_elc;

fftw_complex *h;
fftw_complex *u1_olde,*u2_olde;
fftw_complex *u1_newe,*u2_newe;
fftw_complex *eps_star11e, *eps_star12e, *eps_star22e;
fftw_complex *mu_ele;

fftw_plan planF, planB;

nxny = n_x*n_y;
inv_nxny = 1.0/nxny;

g = fftw_malloc(nxny* sizeof(fftw_complex));
u1_oldc = fftw_malloc(nxny* sizeof(fftw_complex));
u2_oldc = fftw_malloc(nxny* sizeof(fftw_complex));
u1_newc = fftw_malloc(nxny* sizeof(fftw_complex));
u2_newc = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star11c = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star12c = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star22c = fftw_malloc(nxny* sizeof(fftw_complex));
mu_elc = fftw_malloc(nxny* sizeof(fftw_complex));

h = fftw_malloc(nxny* sizeof(fftw_complex));
u1_olde = fftw_malloc(nxny* sizeof(fftw_complex));
u2_olde = fftw_malloc(nxny* sizeof(fftw_complex));
u1_newe = fftw_malloc(nxny* sizeof(fftw_complex));
u2_newe = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star11e = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star12e = fftw_malloc(nxny* sizeof(fftw_complex));
eps_star22e = fftw_malloc(nxny* sizeof(fftw_complex));
mu_ele = fftw_malloc(nxny* sizeof(fftw_complex));

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

//fftw_plan_with_nthreads(8);

planF =
fftw_plan_dft_2d(n_x,n_y,mu_elc,mu_elc,FFTW_FORWARD,FLAG);
planB =
fftw_plan_dft_2d(n_x,n_y,eps_star11c,eps_star11c,FFTW_BACKWARD,FLAG);

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
e = (double *) malloc((size_t) nxny* sizeof(double));
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
J = i2+n_y*i1;
c[J] = creal(comp[J]);
e[J] = creal(eta[J]);
		if(c[J] <= -0.5 || c[J] >= 1.5){
		printf("The composition goes out of bounds1\n");
		printf("5 %d %le\n",J,creal(comp[J]));
		printf("Exiting\n");
		exit(0);
		}

}}

sprintf(NAME,"../output/data/time%dc.dat",
(int) (INDEX + t0));
fpd = fopen(NAME,"w");
tmp = fwrite(&c[0],sizeof(double),(size_t) nxny,fpd);
fclose(fpd);

sprintf(NAME,"../output/data/time%de.dat",
(int) (INDEX + t0));
fpd = fopen(NAME,"w");
tmp = fwrite(&e[0],sizeof(double),(size_t) nxny,fpd);
fclose(fpd);

/** Calculate the Del_sigma_T tensor **/

Del_sigma_Tc = dmatrix(1,2,1,2);
Del_sigma_Te = dmatrix(1,2,1,2);
calculate_Del_sigma_T(DeltaC,epsilon_Tc,Del_sigma_Tc);
calculate_Del_sigma_T(DeltaC,epsilon_Te,Del_sigma_Te);
Ec = dmatrix(1,2,1,2);
strainc = dmatrix(1,2,1,2);
Ee = dmatrix(1,2,1,2);
straine = dmatrix(1,2,1,2);

/* Evolve */


/* Calculate g and its Fourier transform */

fptmp = fopen("../output/noise_info","a");

for(INDEX=0; INDEX < time_steps+1; ++INDEX){

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2+n_y*i1;
realc = creal(comp[J]);
reale = creal(eta[J]);

if(reale < 0){
W = 0.;
Wp = 0.0;
}
else if (reale > 1){
W = 1.0;
Wp = 0.0;
}
else{
W = reale*reale*reale*(10.-15.*reale+6.*reale*reale);
Wp = 3.*reale*reale*(10.-15.*reale+6.*reale*reale)
+ reale*reale*reale*(12.*reale-15.);
}
g[J] = 2.*A*realc + 2.*W*(B*(realc - 1.) - A*realc)+ 0.0*_Complex_I;
h[J] = (B*(1.-realc)*(1.-realc) - A*realc*realc)*Wp +
2.*P*reale*(1.-reale)*(1.-2.*reale);


}}

fftw_execute_dft(planF,g,g);
fftw_execute_dft(planF,h,h);

if(INDEX == 0){

/* The acoustic tensor Omega and the (homogeneous) base solution */

calculate_Omega(n_x,n_y,half_nx,half_ny,kx,ky,Ceff,
Omega11,Omega12,Omega21,Omega22);

calculate_uzero(n_x,n_y,half_nx,half_ny,delta_kx,delta_ky,ave_comp,
comp,Ceff,sigma_Tc,u1_oldc,u2_oldc,n_beta,beta,Omega11,Omega12,Omega21,
Omega22,planF);

calculate_uzero(n_x,n_y,half_nx,half_ny,delta_kx,delta_ky,ave_comp,
eta,Ceff,sigma_Te,u1_olde,u2_olde,n_beta,beta,Omega11,Omega12,Omega21,
Omega22,planF);

}

#ifdef APPROX_CALC
printf("I am doing approximate calculation.\n");
printf("So, I am calculating uzero everytime\n");
if(INDEX != 0){
calculate_uzero(n_x,n_y,half_nx,half_ny,delta_kx,delta_ky,ave_comp,
comp,Ceff,sigma_Tc,u1_oldc,u2_oldc,n_beta,beta,Omega11,Omega12,Omega21,
Omega22,planF);
calculate_uzero(n_x,n_y,half_nx,half_ny,delta_kx,delta_ky,ave_comp,
eta,Ceff,sigma_Te,u1_olde,u2_olde,n_beta,beta,Omega11,Omega12,Omega21,
Omega22,planF);
}
#endif

/* Refine the displacements */

refine_u(n_x,n_y,half_nx,half_ny,delta_x,delta_y,epsilon_Tc,
delta_kx,delta_ky,comp,ave_comp,Ceff,DeltaC,sigma_Tc,
u1_oldc,u2_oldc,u1_newc,u2_newc,Del_sigma_Tc,sig_app,Ec,MAXITR,MAXERR,
n_alpha,n_beta,alpha,beta,Omega11,Omega12,Omega21,Omega22,planF,planB);

refine_u(n_x,n_y,half_nx,half_ny,delta_x,delta_y,epsilon_Te,
delta_kx,delta_ky,eta,ave_comp,Ceff,DeltaC,sigma_Te,
u1_olde,u2_olde,u1_newe,u2_newe,Del_sigma_Te,sig_app,Ee,MAXITR,MAXERR,
n_alpha,n_beta,alpha,beta,Omega11,Omega12,Omega21,Omega22,planF,planB);


for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	eps_star11c[J] = _Complex_I*kx[i1]*u1_newc[J];
	eps_star22c[J] = _Complex_I*ky[i2]*u2_newc[J];
	eps_star12c[J] = 0.5*_Complex_I*(kx[i1]*u2_newc[J]+ky[i2]*u1_newc[J]);
	eps_star11e[J] = _Complex_I*kx[i1]*u1_newe[J];
	eps_star22e[J] = _Complex_I*ky[i2]*u2_newe[J];
	eps_star12e[J] = 0.5*_Complex_I*(kx[i1]*u2_newe[J]+ky[i2]*u1_newe[J]);
}}

/* Get the strain back to the real space */

fftw_execute(planB);
fftw_execute_dft(planB,eps_star12c,eps_star12c);
fftw_execute_dft(planB,eps_star22c,eps_star22c);
fftw_execute_dft(planB,eps_star11e,eps_star11e);
fftw_execute_dft(planB,eps_star12e,eps_star12e);
fftw_execute_dft(planB,eps_star22e,eps_star22e);

/* Calculate mu_el and take it to the Fourier space */

fit=fopen("numerical.dat","w");

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){

	J = i2+n_y*i1;
	
	epsstr11c = creal(eps_star11c[J]*inv_nxny);
	epsstr22c = creal(eps_star22c[J]*inv_nxny);
	epsstr12c = creal(eps_star12c[J]*inv_nxny);
	epsstr11e = creal(eps_star11e[J]*inv_nxny);
	epsstr22e = creal(eps_star22e[J]*inv_nxny);
	epsstr12e = creal(eps_star12e[J]*inv_nxny);
	
	realc = creal(comp[J]);
	reale = creal(eta[J]);
	mu_elc[J] = 0.0;
	mu_ele[J] = 0.0;
	a = alpha[(int) (n_alpha*realc)];
	b = beta[(int) (n_beta*realc)];
	
	ap = alpha_prime[(int) (n_alpha*realc)];
	bp = beta_prime[(int) (n_beta*realc)];
	
	strainc[1][1] = 
	Ec[1][1] + epsstr11c - epsilon_Tc[1][1]*b;
	strainc[1][2] = 
	Ec[1][2] + epsstr12c - epsilon_Tc[1][2]*b;
	strainc[2][1] = 
	Ec[2][1] + epsstr12c - epsilon_Tc[2][1]*b;
	strainc[2][2] = 
	Ec[2][2] + epsstr22c - epsilon_Tc[2][2]*b;

	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
	for(i5=1; i5 < 3; ++i5){
	for(i6=1; i6 < 3; ++i6){
	mu_elc[J] = mu_elc[J]
		+0.5*ap*DeltaC[i3][i4][i5][i6]*strainc[i5][i6]*strainc[i3][i4]
		-bp*epsilon_Tc[i3][i4]*strainc[i5][i6]
			*(Ceff[i3][i4][i5][i6]+DeltaC[i3][i4][i5][i6]*a);
	}}}}
	
	a = alpha[(int) (n_alpha*reale)];
	b = beta[(int) (n_beta*reale)];
	
	ap = alpha_prime[(int) (n_alpha*reale)];
	bp = beta_prime[(int) (n_beta*reale)];
	
	straine[1][1] = 
	Ee[1][1] + epsstr11e - epsilon_Te[1][1]*b;
	straine[1][2] = 
	Ee[1][2] + epsstr12e - epsilon_Te[1][2]*b;
	straine[2][1] = 
	Ee[2][1] + epsstr12e - epsilon_Te[2][1]*b;
	straine[2][2] = 
	Ee[2][2] + epsstr22e - epsilon_Te[2][2]*b;

	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
	for(i5=1; i5 < 3; ++i5){
	for(i6=1; i6 < 3; ++i6){
	mu_ele[J] = mu_ele[J]
		+0.5*ap*DeltaC[i3][i4][i5][i6]*straine[i5][i6]*straine[i3][i4]
		-bp*epsilon_Te[i3][i4]*straine[i5][i6]
			*(Ceff[i3][i4][i5][i6]+DeltaC[i3][i4][i5][i6]*a);
	}}}}

/* The following piece is written by Nitin for comparing with the analytical
solution */


s11=0.0;
s12=0.0;
s21=0.0;
s22=0.0;
for(i11=1;i11<3;i11++){
for(i22=1;i22<3;i22++){
s11=s11+ (Ceff[1][1][i11][i22]+(a*DeltaC[1][1][i11][i22]))*(straine[i11][i22]);
s12=s12+ (Ceff[1][2][i11][i22]+(a*DeltaC[1][2][i11][i22]))*(straine[i11][i22]);
s21=s21+ (Ceff[2][1][i11][i22]+(a*DeltaC[2][1][i11][i22]))*(straine[i11][i22]);
s22=s22+ (Ceff[2][2][i11][i22]+(a*DeltaC[2][2][i11][i22]))*(straine[i11][i22]);
}}
if(i2==n_x/2){
if(i1>=n_y/2){
rr=sqrt((i1-n_x/2)*(i1-n_x/2)+(i2-n_x/2)*(i2-n_x/2));
fprintf(fit,"%lf %lf %lf %lf\n", rr/30.0, s11/4.0, s22/4.0, s21/4.0);
}}}}
fclose(fit);
//exit(0);

/* End of Nitin's code for checking */


fftw_execute(planF);
fftw_execute_dft(planF,comp,comp);
fftw_execute_dft(planF,mu_ele,mu_ele);
fftw_execute_dft(planF,eta,eta);

/* Evolve the composition and eta profiles */


for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){

J = i2+n_y*i1;

k2 = kx[i1]*kx[i1] + ky[i2]*ky[i2];
k4 = k2*k2;

inv_denom = 1.0 + 2.0*M*kappac*k4*delta_t;
inv_denom = 1.0/inv_denom;

comp[J] = inv_denom*( comp[J] - M*k2*delta_t*(g[J]+mu_elc[J]) );

inv_denom = 1.0 + 2.0*L*kappae*k2*delta_t;
inv_denom = 1.0/inv_denom;

eta[J] = inv_denom*( eta[J] - L*delta_t*(h[J]+mu_ele[J]) );

}}

/* Get the composition and eta back to real space */

fftw_execute_dft(planB,comp,comp);
fftw_execute_dft(planB,eta,eta);
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
	J = i2+n_y*i1;
	comp[J] = comp[J]*inv_nxny;
	eta[J] = eta[J]*inv_nxny;
}}

/* Once in a while check that the composition is within bounds,
 * and write the output **/

if(INDEX != 0 && INDEX%STEPS ==0){
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
		J = i2+n_y*i1;
		c[J] = creal(comp[J]);
		e[J] = creal(eta[J]);
		if(c[J] <= -0.5 || c[J] >= 1.5){
		printf("The composition goes out of bounds2\n");
		printf("Exiting\n");

		exit(0);
		}
	}}
	
	sprintf(NAME,"../output/data/time%dc.dat",
		(int) (INDEX + t0));
	fpd = fopen(NAME,"w");
	tmp = fwrite(&c[0],sizeof(double),(size_t) nxny,fpd);
	fclose(fpd);
	sprintf(NAME,"../output/data/time%de.dat",
		(int) (INDEX + t0));
	fpd = fopen(NAME,"w");
	tmp = fwrite(&e[0],sizeof(double),(size_t) nxny,fpd);
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
fftw_free(eps_star11c);
fftw_free(eps_star12c);
fftw_free(eps_star22c);
fftw_free(u1_oldc);
fftw_free(u2_oldc);
fftw_free(u1_newc);
fftw_free(u2_newc);
fftw_free(mu_elc);
fftw_free(h);
fftw_free(eps_star11e);
fftw_free(eps_star12e);
fftw_free(eps_star22e);
fftw_free(u1_olde);
fftw_free(u2_olde);
fftw_free(u1_newe);
fftw_free(u2_newe);
fftw_free(mu_ele);

fftw_destroy_plan(planF);
fftw_destroy_plan(planB);

free(c);
free(e);
free(kx);
free(ky);
free_dmatrix(Del_sigma_Tc,1,2,1,2);
free_dmatrix(Ec,1,2,1,2);
free_dmatrix(strainc,1,2,1,2);
free_dmatrix(Del_sigma_Te,1,2,1,2);
free_dmatrix(Ee,1,2,1,2);
free_dmatrix(straine,1,2,1,2);
free_dvector(Omega11,0,nxny);
free_dvector(Omega12,0,nxny);
free_dvector(Omega21,0,nxny);
free_dvector(Omega22,0,nxny);

}
