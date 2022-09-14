#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void){
FILE *fp;
int i;
double S[201];

fp=fopen("area.dat","r");
for(i=1;i<107;++i){
fscanf(fp,"%le",&S[i]);
}
fclose(fp);

i=3;
while(i<56){
printf("%le\n",(S[i+50]-S[i])/5000.);
i=i+1;
}

return 0;
}
