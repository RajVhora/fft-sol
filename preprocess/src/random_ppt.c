#include<stdlib.h>
#include<math.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_math.h>


void random_ppt(int n_x, int n_y, int R, int N, double comp_precipitate,
double ave_comp, double *comp){

int i1,i2;
int *x,*y;
double X,Y;
double comp_matrix;
double nx_by2,ny_by2;
int nn;
int INDEX;
double tmp1,tmp2,c;

gsl_rng * ran_num;
const gsl_rng_type * Taus;

gsl_rng_env_setup();
Taus = gsl_rng_taus2;
ran_num = gsl_rng_alloc (Taus);

R = R*R;
	
nx_by2 = 0.5*n_x;	
ny_by2 = 0.5*n_y;
	
x = (int *) malloc((size_t) N*sizeof(int));
y = (int *) malloc((size_t) N*sizeof(int));

nn = 0;	
	
for(INDEX=0; INDEX<N; ++INDEX){
	X = (double) gsl_rng_uniform(ran_num);
	Y = (double) gsl_rng_uniform(ran_num);
	X = n_x*X;
	Y = n_y*Y;
	x[INDEX] = (int) X;
	y[INDEX] = (int) Y;

	for(i1=0; i1< n_x; ++i1){
	for(i2=0; i2< n_y; ++i2){
	 tmp1 = 1.0*i1 - 1.0*x[INDEX];
	 if(tmp1 < -nx_by2) tmp1 = tmp1 + 1.0*n_x;
	 if(tmp1 >= nx_by2) tmp1 = tmp1 - 1.0*n_x;
	 tmp1 = tmp1*tmp1/R;
	 tmp2 = 1.0*i2 - 1.0*y[INDEX];
	 if(tmp2 < -ny_by2) tmp2 = tmp2 + 1.0*n_y;
	 if(tmp2 >= ny_by2) tmp2 = tmp2 - 1.0*n_y;
	 tmp2 = tmp2*tmp2/R;
	 c = tmp1 + tmp2;
		if( c <= 1.0 && comp[i2+n_y*i1] != comp_precipitate){
		 comp[i2+n_y*i1] = comp_precipitate;
		 nn = nn + 1;
		}
	}}
}
comp_matrix = (n_x*n_y*ave_comp-comp_precipitate*nn)/(n_x*n_x-nn);
for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_y; ++i2){
 if(comp[i2+n_y*i1] != comp_precipitate){
 comp[i2+n_y*i1] = comp_matrix;}
}}

free(x);
free(y);
gsl_rng_free(ran_num);
}
