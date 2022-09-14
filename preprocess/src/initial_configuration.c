#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "../headers/init_functions.h"

void generate_initial_composition_profile(int n_x, int n_y,
double *comp){

FILE *fp;
int i, j, J;

char initial_profile[7];
char ip1[7] = "spinod";
char ip2[7] = "circul";
char ip3[7] = "random";
char ip4[7] = "ellipt";
char ip5[7] = "thinfi";
char ip6[7] = "twocir";
char junk;
double ave_comp, noise_strength;
int SNA, NOISEADDER;
int R,N,h;
double R1,R2;
double comp_matrix, comp_precipitate;

if( (fp = fopen("input/composition_information","r")) == NULL){
printf("Unable to open input_composition_information. Exiting\n");
exit(0);
}
else{
fp = fopen("input/composition_information","r");
}

fgets(initial_profile,6,fp);

if(strncmp(initial_profile,ip1,5)==0){
fscanf(fp,"%c",&junk);
fscanf(fp,"%le%le",&ave_comp,&noise_strength);
fscanf(fp,"%d%d",&SNA,&NOISEADDER);
spinodal(n_x,n_y,ave_comp,comp);
}
else if(strncmp(initial_profile,ip2,5)==0){
fscanf(fp,"%c",&junk);
fscanf(fp,"%d%le%le%le",
&R,&comp_precipitate,&comp_matrix,&noise_strength);
fscanf(fp,"%d%d",&SNA,&NOISEADDER);
circular(n_x,n_y,R,comp_precipitate,comp_matrix,comp);
}
else if(strncmp(initial_profile,ip3,5)==0){
fscanf(fp,"%c",&junk);
fscanf(fp,"%d%d%le%le%le",
&R,&N,&comp_precipitate,&ave_comp,&noise_strength);
fscanf(fp,"%d%d",&SNA,&NOISEADDER);
random_ppt(n_x,n_y,R,N,comp_precipitate,ave_comp,comp);
}
else if(strncmp(initial_profile,ip4,5)==0){
fscanf(fp,"%c",&junk);
fscanf(fp,"%le%le%le%le%le",
&R1,&R2,&comp_precipitate,&comp_matrix,&noise_strength);
fscanf(fp,"%d%d",&SNA,&NOISEADDER);
elliptic(n_x,n_y,R1,R2,comp_precipitate,comp_matrix,comp);
}
else if(strncmp(initial_profile,ip5,5)==0){
fscanf(fp,"%c",&junk);
fscanf(fp,"%d%le%le%le",
&h,&comp_precipitate,&comp_matrix,&noise_strength);
fscanf(fp,"%d%d",&SNA,&NOISEADDER);
thinfilm(n_x,n_y,h,comp_precipitate,comp_matrix,comp);
}
else if(strncmp(initial_profile,ip6,5)==0){
fscanf(fp,"%c",&junk);
fscanf(fp,"%d%le%le%le",
&R,&comp_precipitate,&comp_matrix,&noise_strength);
fscanf(fp,"%d%d",&SNA,&NOISEADDER);
twocirc(n_x,n_y,R,comp_precipitate,comp_matrix,comp);
}
else{
printf("Unable to determine the initial profile. Exiting\n");
exit(0);
}
fclose(fp);

if(fabs(noise_strength) > 0.5){
printf("Noise strength cannot be greater than 0.5. Exiting\n");
exit(0);
}

if( (fp = fopen("input/noise_strength","w")) == NULL){
 printf("Unable to open input/noise_strength. Exiting\n");
 exit(0);
}
else{
 fp = fopen("input/noise_strength","w");
}
fprintf(fp,"%le\n",noise_strength);
fprintf(fp,"%d\n",SNA);
fprintf(fp,"%d\n",NOISEADDER);
fprintf(fp,"\n");
fprintf(fp,"noise_str, SNA, NOISEADDER\n");
fclose(fp);

if( (fp = fopen("../output/init_profile_info","w")) == NULL){
 printf("Unable to open ../output/init_profile_info. Exiting\n");
 exit(0);
}
else{
 fp = fopen("../output/init_profile_info","w");
}
if(strncmp(initial_profile,ip1,5)==0){
fprintf(fp,"Spinodal: Overall composition of %le with noise\n",
ave_comp);
}
else if(strncmp(initial_profile,ip2,5)==0){
fprintf(fp,"Circular precipitate of radius %d\n",R);
fprintf(fp,"c_p = %le; c_m = %le\n",
comp_precipitate,comp_matrix);
}
else if(strncmp(initial_profile,ip3,5)==0){
fprintf(fp,"%d random precipitates (each of radius %d)\n",
N,R);
fprintf(fp,"c_0 = %le\n",ave_comp);
}
else if(strncmp(initial_profile,ip4,5)==0){
fprintf(fp,"Elliptic precipitate of a=%le and b=%le\n",
R1,R2);
fprintf(fp,"c_p = %le; c_m = %le\n",
comp_precipitate,comp_matrix);
}
else if(strncmp(initial_profile,ip5,5)==0){
fprintf(fp,"Thin film of height %d",h);
fprintf(fp,"c_p = %le; c_m = %le\n",
comp_precipitate,comp_matrix);
}
else if(strncmp(initial_profile,ip6,5)==0){
fprintf(fp,"Two circular precipitates of radius %d\n",R);
fprintf(fp,"c_p = %le; c_m = %le\n",
comp_precipitate,comp_matrix);
}
fclose(fp);
if( (fp = fopen("../output/init_prof","w")) == NULL){
 printf("Unable to open ../output/init_prof. Exiting\n");
 exit(0);
}
else{
 fp = fopen("../output/init_prof","w");
}
for(i=0; i < n_x; ++i){
for(j=0; j < n_y; ++j){
J = j + i*n_y;
fprintf(fp,"%d %d %le\n",i,j,comp[J]);
}}
fclose(fp);

}
