void calculate_sigma_T(double Ceff[3][3][3][3], double eps0_11,
double eps0_22, double eps0_12, double **sigma_T, double **epsilon_T){

int i1,i2,i3,i4;

/* Initialise */

sigma_T[1][1] = 0.0;
sigma_T[1][2] = 0.0;
sigma_T[2][1] = 0.0;
sigma_T[2][2] = 0.0;

/* Get the eigenstrain tensor */

epsilon_T[1][1] = eps0_11;
epsilon_T[2][2] = eps0_22;
epsilon_T[1][2] = eps0_12;
epsilon_T[2][1] = eps0_12;

/* Calculate the eigenstress - Corresponding to Ceff */

for(i1=1; i1<3; ++i1){
for(i2=1; i2<3; ++i2){
for(i3=1; i3<3; ++i3){
for(i4=1; i4<3; ++i4){
sigma_T[i1][i2] = sigma_T[i1][i2] 
	+ Ceff[i1][i2][i3][i4]*epsilon_T[i3][i4];
}}}}

}
