/**************************************************************

This function is based on the following paper:

I. Schmidt and D. Gross

Directional coarsening in Ni-base superalloys: analytical results for
an elasticity-based model

Proc. R. Soc. Lond. A (1999), 455, 3085-3106

The purpose of this function is to predict the direction of elongation
of elliptical precipitates when there is an externally applied
uniaxial stress.

***************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void calculate_aleph(double c11m, double c12m, double c44m, double c11p,
double c12p, double c44p, double **sig_app, double **epsilon_T){

FILE *fp;

/*
Note: Our aleph and w correspond to alpha and little omega of Schmidt 
and Gross. p and m refer to the precipitate and the matrix. A is the
anisotropy parameter and nu is the Poisson's ratio. s is the 
compliance tensor - for the state of plane strain
*/

double Ap,nup,Am,num;
double s11m, s12m, s44m, s11p, s12p, s44p;
double w,p,aleph;

/*
First set of equations in Appendix A
*/

num = 0.5*c12m/(c12m+c44m);
Am = 2.0*c44m/(c11m-c12m);
nup = 0.5*c12p/(c12p+c44p);
Ap = 2.0*c44p/(c11p-c12p);

/*
Equations 5.11
*/

s11m = (3.0+Am-4.0*num)/(8.0*c44m);
s12m = (1.0-Am-4.0*num)/(8.0*c44m);
s44m = (1.0+Am)/(2.0*Am*c44m);

s11p = (3.0+Ap-4.0*nup)/(8.0*c44p);
s12p = (1.0-Ap-4.0*nup)/(8.0*c44p);
s44p = (1.0+Ap)/(2.0*Ap*c44p);

/*
Equations 5.14
*/

if(sig_app[1][1] != 0.0){
p = sig_app[1][1]/(c44m*epsilon_T[1][1]);
}
else p = sig_app[1][1]/(c44m*epsilon_T[1][1]);
w = c44p/c44m;

/*
Equation 5.13
*/

aleph = (1.0 - (s11m - s11p/w)*p)/(1.0 - (s12m - s12p/w)*p);

/*
The following conclusions follow from Figure 3. Let us write the
results of our analysis to an output file. Some of the quantities
printed to the output are calculated using equations 5.18 and 5.20 
*/

if( (fp = fopen("../output/SG_Analysis","w"))==NULL){
printf("Unable to open ../output/SG_Analysis.\n");
printf("Exiting from calculate_aleph.c\n");
exit(0);
}
else{
fp = fopen("../output/SG_Analysis","w");
}

if(fabs(aleph) > 1.0){
fprintf(fp,"Plate perpendicular to the load direction\n");
if(sig_app[1][1] != 0.0){
fprintf(fp,"Namely, the y-direction\n");
}
else
fprintf(fp,"Namely, the x-direction\n");
}
else if(fabs(aleph) < 1.0){ 
fprintf(fp,"Rod parallel to the load direction\n");
if(sig_app[1][1] != 0.0){
fprintf(fp,"Namely, the x-direction\n");
}
else
fprintf(fp,"Namely, the y-direction\n");
}
fprintf(fp,"A_m: %le nu_m: %le\n",Am,num);
fprintf(fp,"s11m: %le s12m: %le s44m: %le\n",s11m,s12m,s44m);
fprintf(fp,"A_p: %le nu_p: %le\n",Ap,nup);
fprintf(fp,"s11p: %le s12p: %le s44p: %le\n",s11p,s12p,s44p);
fprintf(fp,"w: %le p: %le alpha: %le\n",w,p,aleph);
fprintf(fp,"p1_star: 0.0 w1_star: %le\n",((1.0+Ap)/(1.0+Am)));
fprintf(fp,"p2_star: %le  w2_star: %le\n",4.0/(1.0-2.0*num),
((1.0-2.0*nup)/(1.0-2.0*num)));
fclose(fp);

}
