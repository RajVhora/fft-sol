#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

int main(void){
double kappa;
double nu;
double c1, c2, c3, c4;
double mu, mu1, lambda, lambda1;
double A[3][3], B[3];
double Omega11, Omega22, Omega12;
double eps[4];
double a, b, c, d;
double complex x, xbar;
double complex f;
double psidp;
double complex denom;
double tmp;
double epsilon11, epsilon22, epsilon12;
int i,j;
double X,Y;
double sigxx, sigyy,sigxy;
double chi,eta;
double aplusb;
double C1;
double complex C2;
double T1, T2, T3a, T3b;
double T4, T5, T6, T7;
double complex T8, T9;
double complex T10, T11;
double theta;
double complex K1, K2;
double sig11, sig22, sig12;
FILE *fp;
fp=fopen("analytical.dat","w");

a = 25.001; /* Major axis of the elliptic precipitate */
b = 25.0; /* Minor axis of the elliptic precipitate */
mu = 800.0; /*  The shear modulus of the matrix phase */
lambda = 1200.0; /* The other Lame's constant for the matrix phase */
mu1 = 0.5*mu; /* The shear modulus of the precipitate phase */
lambda1 = 0.5*lambda; /*The other Lame's constant for the precipitate phase */
Omega11 = 0.02; /* The 11 component of the eigenstrain */
Omega22 = 0.01; /* The 22 component of the eigenstrain */
Omega12 = 0.005; /* The 12 component of the eigenstrain */
nu = lambda/(2.0*(lambda+mu));

kappa = 3.0-4.0*nu;

c1 = b/((1.0+kappa)*(a+b))*(2.0+kappa+((a-b)/(a+b)));
c2 = b/((1.0+kappa)*(a+b))*(2.0-kappa-((a-b)/(a+b)));
c3 = a/((1.0+kappa)*(a+b))*(2.0-kappa+((a-b)/(a+b)));
c4 = a/((1.0+kappa)*(a+b))*(2.0+kappa-((a-b)/(a+b)));

A[1][1] = 2.0*mu*(c1-1.0) + lambda*(c1+c3-1.0) 
	-2.0*mu1*c1 - lambda1*(c1+c3);
A[1][2] = 2.0*mu*c2 + lambda*(c2+c4-1.0) 
	-2.0*mu1*c2 - lambda1*(c2+c4);
A[2][1] = 2.0*mu*c3 + lambda*(c1+c3-1.0) 
	-2.0*mu1*c3 - lambda1*(c1+c3);
A[2][2] = 2.0*mu*(c4-1.0) + lambda*(c2+c4-1.0) 
	-2.0*mu1*c4 - lambda1*(c2+c4);

B[1] = -2.0*mu1*Omega11-lambda1*(Omega11+Omega22);
B[2] = -2.0*mu1*Omega22-lambda1*(Omega11+Omega22);

eps[1] = (A[2][2]*B[1] - A[1][2]*B[2])/(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
eps[2] = (A[2][1]*B[1] - A[1][1]*B[2])/(A[1][2]*A[2][1]-A[1][1]*A[2][2]);

eps[3] = mu1*Omega12/
(mu-(mu-mu1)*(1.0 - ((4.0*a*b)/((1.0+kappa)*(a+b)*(a+b)))));

epsilon11 = eps[1];
epsilon22 = eps[2];
epsilon12 = eps[3];

printf("Equivalent eigen strain: %le %le %le\n",eps[1],eps[2],eps[3]);

c = sqrt(a*a - b*b);
if(a == b) printf("c is zero. You might get some NaNs\n");
aplusb = a + b;

	for(i=0; i < 750; ++i){
		chi = (double) i*0.01;
	for(j=0; j < 22; ++j){
		eta = j*0.314;
		X = c*cosh(chi)*cos(eta);
		Y = c*sinh(chi)*sin(eta);
		K1 = exp(chi)*(cos(eta)+_Complex_I*sin(eta))
		 - exp(-chi)*(cos(eta)-_Complex_I*sin(eta));
		K2 = exp(chi)*(cos(eta)-_Complex_I*sin(eta))
		 - exp(-chi)*(cos(eta)+_Complex_I*sin(eta));
		theta = 0.5*acos(creal(K1/K2));
		if( (X*X/(a*a) + Y*Y/(b*b)) <= 1.0 ){
			sigxy = 0.0;
			sigyy = (b+2.0*a)*epsilon11/aplusb + b*epsilon22/aplusb;
			sigyy = -4.0*mu*a*sigyy/((1.0+kappa)*aplusb);
			sigxx = a*epsilon11/aplusb + (a+2.0*b)*epsilon22/aplusb;
			sigxx = -4.0*mu*b*sigxx/((1.0+kappa)*aplusb);
			sig22 = 0.5*(sigxx + sigyy - (sigyy-sigxx)*cos(2.0*theta) 
				-2.0*sigxy*sin(2.0*theta));
			sig11 = 0.5*(sigxx + sigyy + (sigyy-sigxx)*cos(2.0*theta) 
				+2.0*sigxy*sin(2.0*theta));
			/*sig12 = 0.5*(2.0*sigxy*cos(2.0*theta) 
				-(sigyy-sigxx)*sin(2.0*theta)); */
			sig12 = -2.0*(mu*a*b*eps[3])/((1.0-nu)*(a+b)*(a+b));
		}
		else{
			C1 = 1.0 - (sinh(2.0*chi)/(cosh(2.0*chi)-cos(2.0*eta)));
			C1 = 8.0*mu*a*b*(epsilon11-epsilon22)*C1/((1.0+kappa)*c*c);
			T1 = cosh(2.0*chi) - ((a*a+b*b)/(c*c));
			T2 = sinh(2.0*chi);
			T3a = 2.0*a*a/(c*c);
			T3b = 2.0*b*b/(c*c);
			T4 = exp(chi);
			T5 = cos(eta);
			T6 = sin(eta);
			T7 = exp(-chi);
			T8 = T4*(T5+_Complex_I*T6) + T7*(T5-_Complex_I*T6);
			T8 = T8/(T4*(T5+_Complex_I*T6) - T7*(T5-_Complex_I*T6));
			T9 = 1.0-exp(-2.0*chi)*cos(2.0*eta) +
			_Complex_I*exp(-2.0*chi)*sin(2.0*eta);
			T10 = T1*T8 - T2 + T3a*T9;
			T11 = -T1*T8 + T2 - T3b*T9;
			C2 = 8.0*mu*a*b/((1.0+kappa)*c*c*(cosh(2.0*chi)-cos(2.0*eta)));
			C2 = C2*(T10*epsilon11+T11*epsilon22);
			sigxy = 0.5*cimag(C2);
			sigxx = 0.5*(C1-creal(C2));
			sigyy = 0.5*(C1+creal(C2));
			sig11 = 0.5*(sigxx + sigyy - (sigyy-sigxx)*cos(2.0*theta) 
				-2.0*sigxy*sin(2.0*theta));
			sig22 = 0.5*(sigxx + sigyy + (sigyy-sigxx)*cos(2.0*theta) 
				+2.0*sigxy*sin(2.0*theta));
			sig12 = 0.5*(2.0*sigxy*cos(2.0*theta) 
				-(sigyy-sigxx)*sin(2.0*theta)); 
		}
		if(j==0){ 
		//printf("%le %le %le %le %le\n",X,Y,sig11,sig22,sig12);
		fprintf(fp,"%le %le %le %le %le\n",X/25.0,Y,sig11/8.0,sig22/8.0,sig12/8.0);
		}
	}}
fclose(fp);
return 0;
}
