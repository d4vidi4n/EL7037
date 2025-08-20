#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/* #include "nrutil.h"*/

void indexx(unsigned long n, double arr[], unsigned long indx[]);

void cont(double *MdistW,double *MdistZ, int k, double *q, unsigned long N){
	double *din,*dout;
	unsigned long *Sdin,*Sdout;
	double aux = 0.0;
	unsigned long i,j,l,h,H,I;
    double A = N*k*(2*N-3*k-1);
    A = 2/A;
	
	din	    = (double *)calloc(N-1,sizeof(double));
	dout	= (double *)calloc(N-1,sizeof(double));
	Sdin	= (unsigned long *)calloc(N-1,sizeof(unsigned long));
	Sdout	= (unsigned long *)calloc(N-1,sizeof(unsigned long));

	for(j=0;j<N;j++){
		l = 0;
		for(i=0;i<N;i++){
			if(i != j){
				din[l]  = *(MdistW+i+N*j);
				dout[l] = *(MdistZ+i+N*j);
				l++;
			}
		}
		indexx(N-1,din,Sdin);
		indexx(N-1,dout,Sdout);

		for(i=1;i<=k;i++){
			I = Sdin[i-1];
			for(h=1;h<=N-1;h++){
				H = Sdout[h-1];
				if(H==I)
					if(h>k)
						aux+= (h-k);
			}
		}
	}
	free(Sdin);
	free(Sdout);
	free(din);
	free(dout);
	*q = 1-A*aux;
}

void qm(double *MdistW,double *MdistZ,int n, int k, double *q, unsigned long N){
	double *din,*dout;
	unsigned long *Sdin,*Sdout;
	double aux = 0.0;
	unsigned long i,j,l,h,H,I;
	
	din	    = (double *)calloc(N-1,sizeof(double));
	dout	= (double *)calloc(N-1,sizeof(double));
	Sdin	= (unsigned long *)calloc(N-1,sizeof(unsigned long));
	Sdout	= (unsigned long *)calloc(N-1,sizeof(unsigned long));

	for(j=0;j<N;j++){
		l = 0;
		for(i=0;i<N;i++){
			if(i != j){
				din[l]  = *(MdistW+i+N*j);
				dout[l] = *(MdistZ+i+N*j);
				l++;
			}
		}
		indexx(N-1,din,Sdin);
		indexx(N-1,dout,Sdout);

		for(i=1;i<=n;i++){
			I = Sdin[i-1];
			for(h=1;h<=N-1;h++){
				H = Sdout[h-1];
				if(H==I)
					if(h==i)
						aux+=3;
					else
						if(h<=n && h>=1 && h!=i)
							aux+=2;
						else
							if(h>n && h<=k && n<k)
								aux+=1;
							else
								aux+=0;
			}
		}
	}
	free(Sdin);
	free(Sdout);
	free(din);
	free(dout);
	*q = aux/(3*N*n);
}

void trust(double *MdistW,double *MdistZ, int k, double *q, unsigned long N){
	double *din,*dout;
	unsigned long *Sdin,*Sdout;
	double aux = 0.0;
	unsigned long i,j,l,h,H,I;
    double A = N*k*(2*N-3*k-1);
    A = 2/A;
	
	din	    = (double *)calloc(N-1,sizeof(double));
	dout	= (double *)calloc(N-1,sizeof(double));
	Sdin	= (unsigned long *)calloc(N-1,sizeof(unsigned long));
	Sdout	= (unsigned long *)calloc(N-1,sizeof(unsigned long));

	for(j=0;j<N;j++){
		l = 0;
		for(i=0;i<N;i++){
			if(i != j){
				din[l]  = *(MdistW+i+N*j);
				dout[l] = *(MdistZ+i+N*j);
				l++;
			}
		}
		indexx(N-1,din,Sdin);
		indexx(N-1,dout,Sdout);

		for(i=1;i<=k;i++){
			I = Sdout[i-1];
			for(h=1;h<=N-1;h++){
				H = Sdin[h-1];
				if(H==I)
					if(h>k)
						aux+= (h-k);
			}
		}
	}
	free(Sdin);
	free(Sdout);
	free(din);
	free(dout);
	*q = 1-A*aux;
}

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(unsigned long n, double arr[], unsigned long indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;
	
	arr--;
	indx--;

	istack=(int *)calloc(NSTACK,sizeof(int));
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) fprintf(stderr,"NSTACK too small in indexx.\n");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free(istack);
	for(j=1;j<=n;j++)
		indx[j]--;
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

