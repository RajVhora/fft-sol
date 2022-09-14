#include<stdio.h>
#include<stdlib.h>

void generate_elasticity_tensor(double c11m, double c12m, double c44m, 
double c11p, double c12p, double c44p, double Ceff[3][3][3][3], 
double DeltaC[3][3][3][3], int alpha_interptype, int prob_type,
double ave_comp){

int i1,i2,i3,i4;

double Cprecip[3][3][3][3];
double Cmatrix[3][3][3][3];

/* Initialise */

for(i1=1; i1 < 3; ++i1){
for(i2=1; i2 < 3; ++i2){
for(i3=1; i3 < 3; ++i3){
for(i4=1; i4 < 3; ++i4){
	Cprecip[i1][i2][i3][i4] = 0.0;
	Cmatrix[i1][i2][i3][i4] = 0.0;
	Ceff[i1][i2][i3][i4] = 0.0;
	DeltaC[i1][i2][i3][i4] = 0.0;
}}}}

if(prob_type == 1){
/* Only problem types 2 (homogeneous elasticity) and 3 (inhomogeneous
 * elasticity) require the calculation of elastic constants */
}
else if(prob_type == 2){
/* There is just a single set of constants */
Ceff[1][1][1][1] = c11m;
Ceff[1][1][2][2] = c12m;
Ceff[1][2][1][2] = c44m;
Ceff[1][2][2][1] = Ceff[2][1][1][2] = Ceff[2][1][2][1]  =
Ceff[1][2][1][2];
Ceff[2][2][1][1] = Ceff[1][1][2][2];
Ceff[2][2][2][2] = Ceff[1][1][1][1];
}
else if(prob_type == 3 || prob_type ==4){
/* Calculate the matrix and precipitate elastic constants */
Cmatrix[1][1][1][1] = c11m;
Cmatrix[1][1][2][2] = c12m;
Cmatrix[1][2][1][2] = c44m;
Cmatrix[1][2][2][1] = Cmatrix[2][1][1][2] = Cmatrix[2][1][2][1]  =
Cmatrix[1][2][1][2];
Cmatrix[2][2][1][1] = Cmatrix[1][1][2][2];
Cmatrix[2][2][2][2] = Cmatrix[1][1][1][1];

Cprecip[1][1][1][1] = c11p;
Cprecip[1][1][2][2] = c12p;
Cprecip[1][2][1][2] = c44p;
Cprecip[1][2][2][1] = Cprecip[2][1][1][2] = Cprecip[2][1][2][1]  =
Cprecip[1][2][1][2];
Cprecip[2][2][1][1] = Cprecip[1][1][2][2];
Cprecip[2][2][2][2] = Cprecip[1][1][1][1];

/* Determine the 'effective' elastic constant Ceff, and the
 * 'perturbative' constant DeltaC - depending on the interpolation
 * funtion used */ 

if(alpha_interptype == 1 || alpha_interptype == 4){
	for(i1=1; i1 < 3; ++i1){
	for(i2=1; i2 < 3; ++i2){
	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
		Ceff[i1][i2][i3][i4] = 
			0.5*(Cprecip[i1][i2][i3][i4] + Cmatrix[i1][i2][i3][i4]);
		DeltaC[i1][i2][i3][i4] = 
			Cprecip[i1][i2][i3][i4] - Cmatrix[i1][i2][i3][i4];
	}}}}
}
else if(alpha_interptype == 2){
	for(i1=1; i1 < 3; ++i1){
	for(i2=1; i2 < 3; ++i2){
	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
		Ceff[i1][i2][i3][i4] = 
			Cmatrix[i1][i2][i3][i4];
		DeltaC[i1][i2][i3][i4] = 
			Cprecip[i1][i2][i3][i4] - Cmatrix[i1][i2][i3][i4];
	}}}}
}
else if(alpha_interptype == 3){
	for(i1=1; i1 < 3; ++i1){
	for(i2=1; i2 < 3; ++i2){
	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
		Ceff[i1][i2][i3][i4] = 
			ave_comp*Cprecip[i1][i2][i3][i4] 
				+ (1.0-ave_comp)*Cmatrix[i1][i2][i3][i4];
		DeltaC[i1][i2][i3][i4] = 
			Cprecip[i1][i2][i3][i4] - Cmatrix[i1][i2][i3][i4];
	}}}}
}
else if(alpha_interptype == 5){
	for(i1=1; i1 < 3; ++i1){
	for(i2=1; i2 < 3; ++i2){
	for(i3=1; i3 < 3; ++i3){
	for(i4=1; i4 < 3; ++i4){
		Ceff[i1][i2][i3][i4] = 
			Cprecip[i1][i2][i3][i4];
		DeltaC[i1][i2][i3][i4] = 
			Cprecip[i1][i2][i3][i4] - Cmatrix[i1][i2][i3][i4];
	}}}}
}
else{
printf("Unable to determine the alpha interpolation type.\n");
printf("Exiting from generate_elast_tensors.c\n");
exit(0);
}
}
else{
printf("Unable to determine the problem type\n");
printf("Exiting from generate_elast_tensors.c\n");
exit(0);
}

}
