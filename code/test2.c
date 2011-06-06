/*
 * test2.c
 *
 */
#include <stdio.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "math.h"
#include "ils.h"
#include "mils.h"
#include <time.h>
#include <stdlib.h>
#include "randn.h"
#include "meschach/matlab.h"
#include "permutationreduction.h"
#include <limits.h>
#include "qrmcp.h"
#include "reduction.h"
#include "search.h"

int main()
{
	MAT *A;
	VEC *z,*y,*v,*low,*up;
	time_t t,start,stop;
	double diff;
	int m,n,i,j;
	double sigma = 0.5;
	m= 10;
	n=7;

	srand(time(NULL));
	y=v_get(m);
	z=v_get(n);
	v=v_get(m);
	low=v_get(n);
	up=v_get(n);
	
	smrand(time(NULL));

	for(i=0;i<n;i++)
	{
		low->ve[i] = INT_MIN;
		up->ve[i] = INT_MAX;
	}

	FILE *fp,*fp2;
	char *matName = "A";
	char *vecName = "y";
	MAT *yMat;
	//If some data exists on disk, load it, otherwise generate some random matrix and vector
	if((fp=fopen("A.mat","rb")) != NULL && (fp2=fopen("y.mat","rb")) != NULL)
	{
		A = m_load(fp,&matName);
		yMat = m_load(fp2,&vecName);
		y = get_col(yMat,0,VNULL);
	}
	else
	{
		A = m_get(m,n);
		m_rand(A);
	

		v_rand(z);
		for(i=0;i<n;i++)
		{
			z->ve[i]=round(10*z->ve[i]);
		}

		for(i=0;i<m;i++)
		{
			v->ve[i]=randn()*sigma;
		}


		v_add(mv_mlt(A,z,VNULL), v, y);
		
		//Try and save input matrix and vector to disk for import to matlab
		if((fp=fopen("A.mat","wb")) != NULL)
		{
			m_save(fp,A,matName);
		}
		if((fp=fopen("y.mat","wb")) != NULL)
		{
			v_save(fp,y,vecName);
		}
	}
	MAT *Z, *R;
	VEC *z_tilt,*x;
	x = v_get(n);
	z_tilt=v_get(n);
	Z = m_get(n, n);
	R = m_get(m,n+1);
	
	_m_copy(A,R,0,0);
	set_col(R,n,y);
	x =v_get(n);
	
	start = time(NULL);
	/*	 Apply QR factorization with minimum column pivoting on A and y to get R.
	 R=[R,y]=[Q^T*A,Q^T*y]*/
	qrmcp(R,Z,n);
	m_resize(R,n,n+1);
   reduction(R, Z, 1, n);
 
   VEC *yp = v_get(m);
	yp = get_col(R,n,VNULL);
	m_resize(R,n,n);

	//Find the permutations, permute R, restore to upper-triangular
	PERM *P = permutationreduction(R,yp,low,up);
	R = px_cols(P,R,MNULL);
	VEC *diag = v_get(n);
	v_output(yp);
	QRfactorY(R,diag,yp);
	//Apply householder transformations to y to get Q'y
	v_output(yp);
	VEC *u = v_get(n);
	double beta,iprod;

	m_resize(R,n,n+1);
	set_col(R,n,yp);
	search(R,n,z_tilt);
	mv_mlt(Z,z_tilt,x);
	stop=time(NULL);
	px_output(P);
	diff=difftime(stop,start);
	printf("Computation time is %e\n",diff);

	M_FREE(A);
	V_FREE(z);
	V_FREE(v);
	V_FREE(y);
	M_FREE(R);
	M_FREE(Z);
	V_FREE(z_tilt);
	V_FREE(x);
	V_FREE(yp);
	return 0;
}
