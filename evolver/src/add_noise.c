#include<stdlib.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_math.h>

void add_noise(int n_x, int n_y, double *comp, double noise_str){

int i1, i2;
double noise, mean_noise;
gsl_rng * ran_num;
const gsl_rng_type * Taus;

gsl_rng_env_setup();
Taus = gsl_rng_taus2;
ran_num = gsl_rng_alloc (Taus);

/* Add a noise of strength noise_str. Keep track of the total 
 * noise added */

mean_noise = 0.0;
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
noise = noise_str*(gsl_rng_uniform(ran_num)-0.5);
mean_noise = mean_noise + noise;
if(noise_str != 0.0){ 
	comp[i2+n_y*i1] = comp[i2+n_y*i1] + noise;}
}}
mean_noise = mean_noise/(n_x*n_y);

/* Make sure that the noise introduction has not changed the overall
 * (average) composition */

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
comp[i2+n_y*i1] = comp[i2+n_y*i1] - mean_noise; 
}}

gsl_rng_free(ran_num);

}
