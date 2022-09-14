#include<stdio.h>
#include<stdlib.h>

#include "gsl/gsl_math.h"

void generate_datafile(char *NAME, char *name, int n_x, int n_y, int N, int st){

FILE *fp;
FILE *FP;
int i,j;
int temp;
int c1,c2,c3,c4,c5;
double r1,r2,r3,r4;
double p1, p2,p3,p4,p5,p6;
double S, S2;
size_t tmp;
double *c;

int half_nx, half_ny;

c=(double *) malloc((size_t)n_x*n_y*sizeof(double));


half_nx = (int) (n_x/2.0);
half_ny = (int) (n_y/2.0);

/** Let us open the data file and read the binary data **/

if( (fp = fopen(NAME,"r")) == NULL){
printf("Unable to open the data file to read.\n");
printf("Exiting from the functions generate_psfile%s\n",name);
exit(0);
}
else{
fp = fopen(NAME,"r");
}
tmp = fread(&c[0],sizeof(double),(size_t) n_x*n_y,fp);
fclose(fp);

/** Let us open the ASCII data file to be written **/

if( (fp = fopen(name,"w")) == NULL){
printf("Unable to open the ps file to write.Exiting\n");
exit(0);
}
else{
fp = fopen(name,"w");
}

if( (FP = fopen("../../../postprocess/data/area.dat","a")) == NULL){
printf("Unable to open the data file to write.Exiting\n");
exit(0);
}
else{
FP = fopen("../../../postprocess/data/area.dat","a");
}



/** Write the ASCII data file **/

c1  = 0;
c2  = 0;
c3  = 0;
c4  = 0;
c5  = 0;
temp = 0;
for(i=0; i < n_x; ++i){
	for(j=0; j < n_y; ++j){
		if(c[j+n_y*i] > 0.5) temp = temp + 1;
		if(i==half_nx){
		 if(c1 == 0 && c[j+n_y*i] < 0.5 && c[j+1+n_y*i] > 0.5){
		 c1 = j;
		 r1 = c[j+n_y*i];
		 c2 = j+1;
		 r2 = c[j+1+n_y*i];
		 c5 = 1;
		 }
		 else if(c5 == 1 && c[j+n_y*i] > 0.5 && c[j+1+n_y*i] < 0.5){
		 c3 = j;
		 r3 = c[j+n_y*i];
		 c4 = j+1;
		 r4 = c[j+1+n_y*i];
		 }
		}
}}
p1 =(0.5-r1)*(c2-c1)/(r2-r1) + c1;
p2 =(0.5-r3)*(c4-c3)/(r4-r3) + c3;

c1  = 0;
c2  = 0;
c3  = 0;
c4  = 0;
c5  = 0;
for(i=0; i < n_x; ++i){
	for(j=0; j < n_y; ++j){

		fprintf(fp,"%d %d %lf\n",i,j,c[j+n_y*i]);

		if(c1 == 0 && c[j+n_y*i] < 0.5 && c[j+n_y*(i+1)] > 0.5){
		 c1 = i;
		 r1 = c[j+n_y*i];
		 c2 = i+1;
		 r2 = c[j+n_y*(i+1)];
		 c5 = 1;
		 }
		 else if(c5 == 1 && c[j+n_y*i] > 0.5 && c[j+n_y*(i+1)] < 0.5){
		 c3 = i;
		 r3 = c[j+n_y*i];
		 c4 = i+1;
		 r4 = c[j+n_y*(i+1)];
		 }

}
/*
fprintf(fp,"\n");
*/
}
p3 =(0.5-r1)*(c2-c1)/(r2-r1) + c1;
p4 =(0.5-r3)*(c4-c3)/(r4-r3) + c3;

for(i=0; i < n_x; ++i){

	for(j=0; j < n_y; ++j){
		if(c[j+n_y*i] > 0.5) temp = temp + 1;
		if(i==j){
		 if(c1 == 0 && c[j+n_y*i] < 0.5 && c[j+1+n_y*i] > 0.5){
		 c1 = j;
		 r1 = c[j+n_y*i];
		 c2 = j+1;
		 r2 = c[j+1+n_y*i];
		 c5 = 1;
		 }
		 else if(c5 == 1 && c[j+n_y*i] > 0.5 && c[j+1+n_y*i] < 0.5){
		 c3 = j;
		 r3 = c[j+n_y*i];
		 c4 = j+1;
		 r4 = c[j+1+n_y*i];
		 }
		}
}}
p5 =(0.5-r1)*(c2-c1)/(r2-r1) + c1;
p6 =(0.5-r3)*(c4-c3)/(r4-r3) + c3;




S = 0.5*((p4-p3)+(p2-p1));
S = 0.5*S;
S2 = S*S;

fprintf(FP,"%d %le %le %le %le %le %le %d\n",st,0.5*(p4-p3),0.5*(p2-p1),S,1.*(p4-p3)/(p2-p1),1.*(p6-p5)/(p2-p1),S2,
tmp);

fclose(fp);
fclose(FP);

free(c);
}

