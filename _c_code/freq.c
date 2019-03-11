#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define m 22.989768*1.6605402e-27
#define mu 9.27400915e-24
#define mum 121.47

int main (int argc, char *argv[]) {
	
	FILE *param, *freq;
	double **par,fl,ft;
	int i,j,N;
	char nome[32],out[32];

	strcpy(nome,"parametri_fit/par_fit_");
	strcat(nome,argv[2]);
	strcat(nome,".dat");	
	param = fopen(nome,"r");

	strcpy(out,"dati/freq_");
	strcat(out,argv[2]);
	strcat(out,".dat");	
	freq = fopen(out,"w+");
	
	N = atoi(argv[1]);
	par = (double **) malloc(N*sizeof(double *));
	for(i=0;i<N;i++) {
		par[i] = (double *) malloc(3*sizeof(double));
	}
	
	for(i=0;i<N;i++) {
		for(j=0;j<3;j++) {
			fscanf(param,"%lf",&par[i][j]);
		}
	}

	for(i=0;i<N;i++) {
		for(j=0;j<3;j++) {
			printf("%lf\n",par[i][j]);
		}
	}

	for(i=0;i<N;i++) {
		fl = sqrt((mum)*par[i][2]);
//		ft = sqrt((mum)*(par[i][1]*par[i][1]/((argv[2][0]=='C'? -1.0 : +1.0 )*par[i][0])-par[i][2]/2.0));
		ft = sqrt((mum)*(par[i][1]*par[i][1]/par[i][0]-par[i][2]/2.0));
		printf("%lf	%15.14lf	%15.14lf\n",10.0+i*10.0,fl,ft);
		fprintf(freq,"%lf	%15.14lf	%15.14lf\n",10.0+i*10.0,fl,ft);
	}
	
	return 0;
}
