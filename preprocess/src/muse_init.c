#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "../headers/init_functions.h"
#include "../headers/nrutil.h"

int main(void){

FILE *fpr, *fpw;

int n_x, n_y;
double delta_x, delta_y;

double kappa1,kappa2;
double M,L;
double A,B,P;

double alpha,beta,gamma,pre_factor;

int time_steps1, time_steps2, time_steps3;
double delta_t1, delta_t2, delta_t3;
double t0;

double *comp, noise_str;
int SNA, NOISEADDER;

int STEPS;

/* Read the data from the files in the input directory */

if( (fpr = fopen("input/system_size","r")) == NULL){
printf("Unable to open input/system_size. Exiting");
exit(0);
}
else{
fpr = fopen("input/system_size","r");
}
fscanf(fpr,"%d%d%le%le",&n_x,&n_y,&delta_x,&delta_y);
fclose(fpr);

if( (fpr = fopen("input/constants","r")) == NULL){
printf("Unable to oepn input/constants. Exiting");
exit(0);
}
else{
fpr = fopen("input/constants","r");
}
fscanf(fpr,"%le%le",&kappa1,&kappa2);
fscanf(fpr,"%le%le",&M,&L);
fscanf(fpr,"%le%le%le\n",&A,&B,&P);
fclose(fpr);

if( (fpr = fopen("input/time_information","r")) == NULL){
printf("Unable to open input/time_information. Exiting");
exit(0);
}
else{
fpr = fopen("input/time_information","r");
}
fscanf(fpr,"%le",&t0);
fscanf(fpr,"%le%d",&delta_t1,&time_steps1);
fscanf(fpr,"%le%d",&delta_t2,&time_steps2);
fscanf(fpr,"%le%d",&delta_t3,&time_steps3);
fclose(fpr);

if( (fpr = fopen("input/magic_number","r")) == NULL){
printf("Unable to open input/magic_number. Exiting");
exit(0);
}
else{
fpr = fopen("input/magic_number","r");
}
fscanf(fpr,"%d",&STEPS);
fclose(fpr);

if( (fpr = fopen("input/noise_strength","r")) == NULL){
printf("Unable to open input/noise_strength. Exiting");
exit(0);
}
else{
fpr = fopen("input/noise_strength","r");
}
fscanf(fpr,"%le",&noise_str);
fclose(fpr);

/**
Write the data to a file
**/

if( (fpw = fopen("../input/system_constants","w")) == NULL){
printf("Unable to open ../input/system_constants. Exiting");
exit(0);
}
else{
fpw = fopen("../input/system_constants","w");
}
fprintf(fpw,"%d\n%d\n",n_x,n_y);
fprintf(fpw,"%le\n%le\n",delta_x,delta_y);
fprintf(fpw,"%le\n%le\n",kappa1,kappa2);
fprintf(fpw,"%le\n%le\n%le\n",A,B,P);
fprintf(fpw,"%le\n",t0);
fprintf(fpw,"%d\t%le\n",time_steps1,delta_t1);
fprintf(fpw,"%d\t%le\n",time_steps2,delta_t2);
fprintf(fpw,"%d\t%le\n",time_steps3,delta_t3);
fprintf(fpw,"%d\n\n",STEPS);
fprintf(fpw,"\n\nn_x,n_y,delta_x,delta_y,kappac,kappae,A\n");
fprintf(fpw,"B,P\n");
fprintf(fpw,"time_steps,delta_t,STEPS\n");
fclose(fpw);

/* Generate the initial composition profile */

comp = dvector(0,n_x*n_y-1);
generate_initial_composition_profile(n_x,n_y,comp);

if( (fpr = fopen("input/noise_strength","r")) == NULL){
printf("Unable to open input/noise_strength. Exiting");
exit(0);
}
else{
fpr = fopen("input/noise_strength","r");
}
fscanf(fpr,"%le",&noise_str);
fscanf(fpr,"%d",&SNA);
fscanf(fpr,"%d",&NOISEADDER);
fclose(fpr);

/* Add noise to the initial profile */

add_noise(n_x,n_y,comp,noise_str,SNA,NOISEADDER);  

free_dvector(comp,0,n_x*n_y-1);

return 0;
}
