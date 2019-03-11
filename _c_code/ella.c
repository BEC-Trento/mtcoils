#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include<gsl/gsl_sf_ellint.h>

#include "nrutil.h"//copiare il file dal ricettario
#include "nrutil.c"
#define ERRTOL_f 0.00025 //parametri per le funzioni ellittiche
#define TINY_f 1.5e-38
#define BIG_f 3.0e37
#define THIRD (1.0/3.0)
#define C1_f (1.0/24.0)
#define C2_f 0.1
#define C3_f (3.0/44.0)
#define C4_f (1.0/14.0)
#define ERRTOL_d 0.00015
#define TINY_d 1.0e-25
#define BIG_d 4.5e21
#define C1_d (3.0/14.0)
#define C2_d (1.0/6.0)
#define C3_d (9.0/22.0)
#define C4_d (3.0/26.0)
#define C5_d (0.25*C3_d)
#define C6_d (1.5*C4_d)

#define CON 1.4 //parametri per la derivazione numerica
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

#define pi 3.141592654
#define mu 1.257e-6

#define N 1000 //punti della grigla di calcolo
#define inf -0.005 //limiti della regione di calcolo
#define sup 0.005

//per introdurre i parametri di bobina userei una struct
	typedef struct {
		int spire_l;
		int spire_r;
		double dim_l;
		double dim_r;
		double corrente;
		double dist;
		double raggio_int;
	} bobina;

short int segno (double var) {

	if (var > 0.0)
		return 1;
	else if (var < 0.0)
		return -1;
	else if (var == 0.0)
		return 0;

}

//integrali ellittici
double rf (double x, double y, double z) {

	double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;

	if (DMIN(DMIN(x,y),z) < 0.0 || DMIN(DMIN(x+y,x+z),y+z) < TINY_f || DMAX(DMAX(x,y),z) > BIG_f)
		nrerror("Argomenti non validi in rf\n");

	xt = x;
	yt = y;
	zt = z;
	do {
		sqrtx = sqrt(xt);
		sqrty = sqrt(yt);
		sqrtz = sqrt(zt);
		alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt = 0.25*(xt+alamb);
		yt = 0.25*(yt+alamb);
		zt = 0.25*(zt+alamb);
		ave = THIRD*(xt+yt+zt);
		delx = (ave-xt)/ave;
		dely = (ave-yt)/ave;
		delz = (ave-zt)/ave;
	} while (DMAX(DMAX(abs(delx),abs(dely)),abs(delz)) > ERRTOL_f);
	e2 = delx*dely-delz*delz;
	e3 = delx*dely*delz;

	return (1.0+(C1_f*e2-C2_f-C3_f*e3)*e2+C4_f*e3)/sqrt(ave);
}

double rd (double x, double y, double z) {

	double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;

	if (DMIN(x,y) < 0.0 || DMIN(x+y,z) < TINY_d || DMAX(DMAX(x,y),z) > BIG_d)
		nrerror("Argomenti non validi in rd");

	xt = x;
	yt = y;
	zt = z;
	sum = 0.0;
	fac = 1.0;
	do {
		sqrtx = sqrt(x);
		sqrty = sqrt(y);
		sqrtz = sqrt(z);
		alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac = 0.25*fac;
		xt = 0.25*(xt+alamb);
		yt = 0.25*(yt+alamb);
		zt = 0.25*(zt+alamb);
		ave = 0.2*(xt+yt+3.0*zt);
		delx = (ave-xt)/ave;
		dely = (ave-yt)/ave;
		delz = (ave-zt)/ave;
	} while (DMAX(DMAX(abs(delx),abs(dely)),abs(delz)) > ERRTOL_d);
	ea = delx*dely;
	eb = delz*delz;
	ec = ea-eb;
	ed = ea-6.0*eb;
	ee = ed+ec+ec;

	return 3.0*sum+fac*(1.0+ed*(-C1_d+C5_d*ed-C6_d*delz*ee)+delz*(C2_d*ee+delz*(-C3_d*ec+delz*C4_d*ea)))/(ave*sqrt(ave));
}

double ellf (double phi, double ak) {//integrale ellittico del primo tipo
// F(phi,k), per calcolare K(k)=F(pi/2,k)
	double rf(double x, double y, double z);
	double s;

	s = sin(phi);
	return s*rf(DSQR(cos(phi)),(1.0-s*ak)*(1.0+s*ak),1.0);
}

double elle (double phi, double ak) {//integrale ellittico del secondo tipo
// E(phi,k), per calcolare E(k)=E(pi/2,k)
	double rd(double x, double y, double z);
	double rf(double x, double y, double z);
	double cc,q,s;

	s = sin(phi);
	cc = DSQR(cos(phi));
	q = (1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-(DSQR(s*ak))*rd(cc,q,1.0)/3.0);
}

double Basse (double x, double r, bobina param) {//calcolo componente assiale
//controllare l'allocazione dinamica di b[][]
	int l,m,i;
	double B, **b;
	double a_m, x_l, alpha_m, beta_lm, gamma_l, Q_lm, k_lm;

	b = (double **) malloc(param.spire_l*sizeof(double *));
	for(i=0;i<param.spire_l;i++) {
		b[i] = (double *) malloc(param.spire_r*sizeof(double));
	}

	B = 0.0;
	for(l=0;l<param.spire_l;l++) {
		for(m=0;m<param.spire_r;m++) {
			a_m = param.raggio_int + param.dim_r*(m+0.5);
			x_l = param.dist + segno(param.dist)*param.dim_l*(l+0.5);
			alpha_m = r/a_m;
			beta_lm = (x-x_l)/a_m;
			gamma_l = (x-x_l)/r;
			Q_lm = pow(1+alpha_m,2)+pow(beta_lm,2);
			k_lm = sqrt(4*alpha_m/Q_lm);

			//b[l][m] = mu*param.corrente*0.5/(a_m*pi*sqrt(Q_lm))*(elle(pi/2,k_lm)*(1-pow(alpha_m,2)-pow(beta_lm,2))/(Q_lm-4*alpha_m)+ellf(pi/2,k_lm));
			b[l][m] = mu*param.corrente*0.5/(a_m*pi*sqrt(Q_lm))*(gsl_sf_ellint_Ecomp(k_lm,0)*(1-pow(alpha_m,2)-pow(beta_lm,2))/(Q_lm-4*alpha_m)+gsl_sf_ellint_Kcomp(k_lm,0));
			B = B + b[l][m];
		}
	}

	//disalloca la memoria
	for(i=0;i<param.spire_l;i++) {
		free(b[i]);
	}
	free(b);

	return B;
}

double Braggio (double x, double r, bobina param) {//calcolo componente radiale

	int l,m,i;
	double B, **b;
	double a_m, x_l, alpha_m, beta_lm, gamma_l, Q_lm, k_lm;

	b = (double **) malloc(param.spire_l*sizeof(double *));
	for(i=0;i<param.spire_l;i++) {
		b[i] = (double *) malloc(param.spire_r*sizeof(double));
	}

	B = 0.0;
	for(l=0;l<param.spire_l;l++) {
		for(m=0;m<param.spire_r;m++) {
			a_m = param.raggio_int + param.dim_r*(m+0.5);
			x_l = param.dist + segno(param.dist)*param.dim_l*(l+0.5);
			alpha_m = r/a_m;
			beta_lm = (x-x_l)/a_m;
			gamma_l = (x-x_l)/r;
			Q_lm = pow(1+alpha_m,2)+pow(beta_lm,2);
			k_lm = sqrt(4*alpha_m/Q_lm);

			//b[l][m] = mu*param.corrente*0.5*gamma_l/(a_m*pi*sqrt(Q_lm))*(elle(pi/2,k_lm)*(1+pow(alpha_m,2)+pow(beta_lm,2))/(Q_lm-4*alpha_m)-ellf(pi/2,k_lm));
			b[l][m] = mu*param.corrente*0.5*gamma_l/(a_m*pi*sqrt(Q_lm))*(gsl_sf_ellint_Ecomp(k_lm,0)*(1+pow(alpha_m,2)+pow(beta_lm,2))/(Q_lm-4*alpha_m)-gsl_sf_ellint_Kcomp(k_lm,0));
			B = B + b[l][m];
		}
	}

	//disalloca la memoria
	for(i=0;i<param.spire_l;i++) {
		free(b[i]);
	}
	free(b);

	return B;
}

//derivate numeriche
void dfridr (double *vec, double *dvec, double *d2vec, double h) {

	int i;

	if (h == 0.0) nrerror("h deve essere diverso da zero!\n");
//valori iniziali
	dvec[0] = -(vec[0]-vec[1])/h;
	d2vec[0] = (vec[2]-2.0*vec[1]+vec[0])/(h*h);

//valori intermedi
	for(i=1;i<(N-1);i++) {
		dvec[i] = (vec[i+1]-vec[i-1])/(2.0*h);
		d2vec[i] = (vec[i-1] - 2.0*vec[i] + vec[i+1])/(h*h);
 	}

//valori finali
	dvec[N-1] = (vec[N-1]-vec[N-2])/h;
	d2vec[N-1] = (dvec[N-1]-dvec[N-2])/h;//inserire l'espressione corretta!

}



int main (int argc, char *argv[]) {

	int i,j;
	double pos[N], campi[2][N], grad[2][N], curv[2][N];//0-x e 1-z
	double h;
	double compensazione[2][N];
	char nome[32];
	int flag[2];
	int p[2];

	FILE *totale,*param, *err;

	strcpy(nome,"dati/totale");
	strcat(nome,argv[2]);
	strcat(nome,argv[1]);
	strcat(nome,".dat");
	totale = fopen(nome,"w+");

	strcpy(nome,"dati/parametri");
	strcat(nome,argv[2]);
	strcat(nome,argv[1]);
	strcat(nome,".plt");
	param = fopen(nome,"w+");

	strcpy(nome,"dati/err");
	strcat(nome,argv[2]);
	strcat(nome,argv[1]);
	strcat(nome,".dat");
	err = fopen(nome,"w+");

//bobina superiore quadrupolo
	bobina quadrupolo_h;

	quadrupolo_h.spire_l = 8;
	quadrupolo_h.spire_r = 9;
	quadrupolo_h.dim_l = 0.0042;
	quadrupolo_h.dim_r = 0.0042;
	quadrupolo_h.corrente = -atof(argv[1]);//-200.0;
	quadrupolo_h.dist = 0.02;
	quadrupolo_h.raggio_int = 0.055;

//bobina inferiore quadrupolo
	bobina quadrupolo_l;

	quadrupolo_l.spire_l = 8;
	quadrupolo_l.spire_r = 9;
	quadrupolo_l.dim_l = 0.0042;
	quadrupolo_l.dim_r = 0.0042;
	quadrupolo_l.corrente = atof(argv[1]);//200.0;
	quadrupolo_l.dist = -0.02;
	quadrupolo_l.raggio_int = 0.055;

//bobina curvatura
	bobina pinch;

	pinch.spire_l = 3;
	pinch.spire_r = 4;
	pinch.dim_l = 0.0032;
	pinch.dim_r = 0.0032;
	pinch.corrente = atof(argv[1]);//200.0;
	pinch.dist = -0.0195;
	pinch.raggio_int = 0.0125;

//bobina compensazione
	bobina comp1;

	comp1.spire_l = 4;
	comp1.spire_r = 4;
	comp1.dim_l = 0.0042;
	comp1.dim_r = 0.0042;
	comp1.corrente = -atof(argv[1]);//-200.0;
	comp1.dist = -0.075;
	comp1.raggio_int = 0.087;

	bobina comp2;

	comp2.spire_l = 4;
	comp2.spire_r = 4;
	comp2.dim_l = 0.0042;
	comp2.dim_r = 0.0042;
	comp2.corrente = -atof(argv[1]);//200.0;
	comp2.dist = 0.075;
	comp2.raggio_int = 0.087;

//definisco la griglia sugli assi x (Quadrupolo) e z (Curvatura)
	for(i=0;i<N;i++) {
		pos[i] = inf + i*(sup-inf)/(N-1);
		campi[0][i]=0.0;
		campi[1][i]=0.0;
	}

//calcolo una bobina di quadrupolo
	h = quadrupolo_l.dist + segno(quadrupolo_l.dist)*(quadrupolo_l.spire_l*quadrupolo_l.dim_l)/2;
	for(i=0;i<N;i++) {
		campi[0][i] = ( pos[i] < h ?
					Basse(pos[i],0.0,quadrupolo_l) :
					Basse(pos[i],0.0,quadrupolo_l));
		campi[1][i] = ( pos[i]<0.0 ?
					-Braggio(0.0,-pos[i],quadrupolo_l) :
					Braggio(0.0,pos[i],quadrupolo_l));
	}

//calcolo altra bobina di quadrupolo
	h = quadrupolo_h.dist + segno(quadrupolo_h.dist)*(quadrupolo_h.spire_l*quadrupolo_h.dim_l)/2;
	for(i=0;i<N;i++) {
		campi[0][i] = campi[0][i] + ( pos[i] < h ?
					Basse(pos[i],0.0,quadrupolo_h) :
					Basse(pos[i],0.0,quadrupolo_h));
		campi[1][i] = campi[1][i] + ( pos[i]<0.0 ?
					-Braggio(0.0,-pos[i],quadrupolo_h) :
					Braggio(0.0,pos[i],quadrupolo_h));
	}

//calcolo bobina di curvatura
	h = pinch.dist + segno(pinch.dist)*(pinch.spire_l*pinch.dim_l)/2;
	for(i=0;i<N;i++) {
		campi[0][i] = campi[0][i] + ( pos[i] < 0.0 ?
					-Braggio(0.0,-pos[i],pinch) :
					Braggio(0.0,pos[i],pinch));
		campi[1][i] = campi[1][i] + ( pos[i] < h ?
					Basse(pos[i],0.0,pinch) :
					Basse(pos[i],0.0,pinch));
	}

	if(argv[2][0]=='C') {
	//calcolo bobine di compensazione
		h = comp1.dist + segno(comp1.dist)*(comp1.spire_l*comp1.dim_l)/2;
		for(i=0;i<N;i++) {
			compensazione[0][i] = ( pos[i] < 0.0 ?
						-Braggio(0.0,-pos[i],comp1) :
						Braggio(0.0,pos[i],comp1));
			compensazione[1][i] =  ( pos[i] < h ?
						Basse(pos[i],0.0,comp1) :
						Basse(pos[i],0.0,comp1));
		}
		h = comp2.dist + segno(comp2.dist)*(comp2.spire_l*comp2.dim_l)/2;
		for(i=0;i<N;i++) {
			compensazione[0][i] = compensazione[0][i] + ( pos[i] < 0.0 ?//
							-Braggio(0.0,-pos[i],comp2) ://
							Braggio(0.0,pos[i],comp2));
			compensazione[1][i] = compensazione[1][i] + ( pos[i] < h ?//
						Basse(pos[i],0.0,comp2) ://
						Basse(pos[i],0.0,comp2));
		}

	//sommo la compensazione
		for (i=0;i<2;i++) {
			for (j=0;j<N;j++) {
				campi[i][j] = campi[i][j] + compensazione[i][j];
//				if(i==0)
//				printf("%d # %d ## %15.14lf\n",i,j,compensazione[i][j]);
			}
		}
	}

//calcolo gradiente e curvatura
	h = (sup-inf)/(N-1);
	for (i=0;i<2;i++) {
		dfridr(campi[i],grad[i],curv[i],h);
	}

//
	for (i=0;i<2;i++) {
		for (j=0;j<N;j++) {
			if (segno(grad[i][j])!=segno(grad[i][j-1])) {
				flag[i]=j;
				printf("%d	%d	%lf\n",i,flag[i],pos[flag[i]]);
			}
		}
	}

//errore
	for (i=0;i<(N-1);i++) {
		fprintf(err,"%9.8lf	%5.4lf\n",pos[i],(0.0035+ ((campi[1][flag[1]]+grad[1][N/2-1]*pos[i]+0.5*curv[1][flag[1]]*pos[i]*pos[i])-(campi[1][i]))/campi[1][i] ) );
	}
	fclose(err);


//generazione file campi
	for(i=0;i<N;i++) {
		fprintf(totale,"%9.8lf	%9.8lf	%9.8lf	%9.8lf	%9.8lf	%9.8lf	%9.8lf\n",pos[i],campi[0][i],campi[1][i],grad[0][i],grad[1][i],curv[0][i],curv[1][i]);
	}

	fclose(totale);


//generazione file parametri
	fprintf(param,"b0 = %9.8lf\n",campi[1][flag[1]]);
	fprintf(param,"b1 = %9.8lf\n",grad[1][N/2-1]);
	fprintf(param,"bd = %9.8lf\n",grad[0][N/2]);
	fprintf(param,"b2d = %9.8lf\n",curv[1][flag[1]]);
	fclose(param);

	return 0;
}
