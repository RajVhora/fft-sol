#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <fftw3.h>

#include "../headers/functions.h"

void write_ppfiles(int n_x, int n_y, double delta_t1, int time_steps1, 
double t0, double delta_t2, int time_steps2, double t1, 
double delta_t3, int time_steps3, double t2, int STEPS){

FILE *fpp;

int i1,st;

char NAME[50];
char name[50];

/* Write the psmaker file for c variable */

fpp =  fopen("../output/postprocessor/ps/src/psmakerc.c","a");
fprintf(fpp,"/** This is a machine generated C-file. **/\n");
fprintf(fpp,"#include<stdio.h>\n");
fprintf(fpp,"#include<stdlib.h>\n");
fprintf(fpp,"\n");
fprintf(fpp,"extern void generate_psfile(char *NAME,char *name,\n");
fprintf(fpp,"int n_x,int n_y,int N,int st);\n");
fprintf(fpp,"\n");
fprintf(fpp,"int main(void){\n");
fprintf(fpp,"FILE *fp;\n");
for(i1=0; i1 < time_steps1+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%dc.dat",
(int) (i1 + t0));
sprintf(name,"../../../postprocess/ps/ctime%d.ps", 
(int) (i1 + t0));
st=i1+t0;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the psfile %s\\n\");\n",name);
fprintf(fpp,"generate_psfile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps2+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%dc.dat",
(int) (i1 + t1));
sprintf(name,"../../../postprocess/ps/ctime%d.ps", 
(int) (i1 + t1));
st=i1+t1;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the psfile %s\\n\");\n",name);
fprintf(fpp,"generate_psfile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps3+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%dc.dat",
(int) (i1 + t2));
sprintf(name,"../../../postprocess/ps/ctime%d.ps", 
(int) (i1 + t2));
st=i1+t2;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the psfile %s\\n\");\n",name);
fprintf(fpp,"generate_psfile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
fprintf(fpp,"return(0);\n");
fprintf(fpp,"}\n");
fclose(fpp);

/* Write the datamaker file for c variable*/

fpp =  fopen("../output/postprocessor/data/src/datamakerc.c","a");
fprintf(fpp,"/** This is a machine generated C-file. **/\n");
fprintf(fpp,"#include<stdio.h>\n");
fprintf(fpp,"#include<stdlib.h>\n");
fprintf(fpp,"\n");
fprintf(fpp,"extern void generate_datafile(char *NAME,char *name,\n");
fprintf(fpp,"int n_x,int n_y,int N, int st);\n");
fprintf(fpp,"\n");
fprintf(fpp,"int main(void){\n");
fprintf(fpp,"FILE *fp;\n");
for(i1=0; i1 < time_steps1+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%dc.dat", 
(int) (i1 +  t0));
sprintf(name,"../../../postprocess/data/ctime%d.data", 
(int) (i1 + t0));
st=i1+t0;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the datafile %s\\n\");\n",name);
fprintf(fpp,"generate_datafile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps2+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%dc.dat", 
(int) (i1 +  t1));
sprintf(name,"../../../postprocess/data/ctime%d.data", 
(int) (i1 + t1));
st=i1+t1;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the datafile %s\\n\");\n",name);
fprintf(fpp,"generate_datafile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps3+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%dc.dat", 
(int) (i1 +  t2));
sprintf(name,"../../../postprocess/data/ctime%d.data", 
(int) (i1 + t2));
st=i1+t2;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the datafile %s\\n\");\n",name);
fprintf(fpp,"generate_datafile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
fprintf(fpp,"return(0);\n");
fprintf(fpp,"}\n");
fclose(fpp);

//

/* Write the psmaker file for e variable */

fpp =  fopen("../output/postprocessor/ps/src/psmakere.c","a");
fprintf(fpp,"/** This is a machine generated C-file. **/\n");
fprintf(fpp,"#include<stdio.h>\n");
fprintf(fpp,"#include<stdlib.h>\n");
fprintf(fpp,"\n");
fprintf(fpp,"extern void generate_psfile(char *NAME,char *name,\n");
fprintf(fpp,"int n_x,int n_y,int N, int st);\n");
fprintf(fpp,"\n");
fprintf(fpp,"int main(void){\n");
fprintf(fpp,"FILE *fp;\n");
for(i1=0; i1 < time_steps1+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%de.dat",
(int) (i1 + t0));

sprintf(name,"../../../postprocess/ps/etime%d.ps", 
(int) (i1 + t0));
st=i1+t0;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the psfile %s\\n\");\n",name);
fprintf(fpp,"generate_psfile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps2+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%de.dat",
(int) (i1 + t1));
sprintf(name,"../../../postprocess/ps/etime%d.ps", 
(int) (i1 + t1));
st=i1+t1;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the psfile %s\\n\");\n",name);
fprintf(fpp,"generate_psfile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps3+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%de.dat",
(int) (i1 + t2));
sprintf(name,"../../../postprocess/ps/etime%d.ps", 
(int) (i1 + t2));
st=i1+t2;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the psfile %s\\n\");\n",name);
fprintf(fpp,"generate_psfile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
fprintf(fpp,"return(0);\n");
fprintf(fpp,"}\n");
fclose(fpp);

/* Write the datamaker file for e variable*/

fpp =  fopen("../output/postprocessor/data/src/datamakere.c","a");
fprintf(fpp,"/** This is a machine generated C-file. **/\n");
fprintf(fpp,"#include<stdio.h>\n");
fprintf(fpp,"#include<stdlib.h>\n");
fprintf(fpp,"\n");
fprintf(fpp,"extern void generate_datafile(char *NAME,char *name,\n");
fprintf(fpp,"int n_x,int n_y,int N, int st);\n");
fprintf(fpp,"\n");
fprintf(fpp,"int main(void){\n");
fprintf(fpp,"FILE *fp;\n");
for(i1=0; i1 < time_steps1+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%de.dat", 
(int) (i1 +  t0));
sprintf(name,"../../../postprocess/data/etime%d.data", 
(int) (i1 + t0));
st = i1+t0;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the datafile %s\\n\");\n",name);
fprintf(fpp,"generate_datafile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps2+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%de.dat", 
(int) (i1 +  t1));
sprintf(name,"../../../postprocess/data/etime%d.data", 
(int) (i1 + t1));
st=i1+t1;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the datafile %s\\n\");\n",name);
fprintf(fpp,"generate_datafile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
for(i1=0; i1 < time_steps3+1; i1=i1+STEPS){
sprintf(NAME,"../../data/time%de.dat", 
(int) (i1 +  t2));
sprintf(name,"../../../postprocess/data/etime%d.data", 
(int) (i1 + t2));
st = i1+t2;
fprintf(fpp,"if( (fp=fopen(\"%s\",\"r\")) == NULL){\n",name);
fprintf(fpp,"printf(\"Generating the datafile %s\\n\");\n",name);
fprintf(fpp,"generate_datafile(\"%s\",\"%s\",%d,%d,%d,%d);\n",NAME,name,n_x,n_y,n_x,st);
fprintf(fpp,"}\n");
}
fprintf(fpp,"return(0);\n");
fprintf(fpp,"}\n");
fclose(fpp);

}
