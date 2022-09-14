#include<stdio.h>
#include<stdlib.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_math.h>


void add_noise(int n_x, int n_y, double *comp, 
double noise_str, int SNA, int NOISEADDER){

FILE *fp;

int i1, i2;
double average_composition;
gsl_rng * ran_num;
const gsl_rng_type * Taus;
size_t tmp;

double noise, mean_noise;

double *temp;
int J;

gsl_rng_env_setup();
Taus = gsl_rng_taus2;
ran_num = gsl_rng_alloc (Taus);

temp = (double *) malloc((size_t) n_x*n_y* sizeof(double));
mean_noise = 0.0;
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
noise = noise_str*(gsl_rng_uniform(ran_num)-0.5);
mean_noise = mean_noise + noise;
if(noise_str != 0.0) 
comp[i2+n_y*i1] = comp[i2+n_y*i1] + noise;
}}
mean_noise = mean_noise/(n_x*n_y);

average_composition = 0.0;
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
comp[i2+n_y*i1] = comp[i2+n_y*i1] - mean_noise; 
average_composition = average_composition + comp[i2+n_y*i1];
}}
average_composition = average_composition/(n_x*n_y);

fp = fopen("../input/ave_comp","w");
tmp = fwrite(&average_composition,sizeof(double),(size_t) 1,fp);
fclose(fp);
fp = fopen("../input/init_compo_prof","w");
tmp = fwrite(&comp[0],sizeof(double),(size_t)n_x*n_y,fp);
fclose(fp);
fp = fopen("../input/init_eta_prof","w");
tmp = fwrite(&comp[0],sizeof(double),(size_t)n_x*n_y,fp);
fclose(fp);

/*
fp = fopen("../input/init_eta1_prof","w");
tmp = fwrite(&comp[0],sizeof(double),(size_t)n_x*n_y,fp);
fclose(fp);
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
J = i2 + n_y*i1;
if(comp[J] == 0)
temp[J] = 1.0;
else
temp[J] = 0.0;
}}
fp = fopen("../input/init_eta2_prof","w");
tmp = fwrite(&temp[0],sizeof(double),(size_t)n_x*n_y,fp);
fclose(fp);

fp = fopen("../input/sustained_noise","w");
fprintf(fp,"%d\n",SNA);
fprintf(fp,"%d\n",NOISEADDER);
fprintf(fp,"%le\n",noise_str);
fclose(fp);
fp = fopen("../output/randomnumber_info","w");
fprintf(fp,"Random number generator: %s\n",gsl_rng_name (ran_num));
fprintf(fp,"Seed: %lu\n",gsl_rng_default_seed);
fclose(fp);
*/


gsl_rng_free(ran_num);
free(temp);
}
