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

int main()
{
	MAT *A;
	VEC *z,*y,*v,*l,*u;
	time_t t,start,stop;
	double diff;
	int m,n,i;
	double sigma = 0.5;
	m=10;
	n=7;

	srand(time(NULL));
	y=v_get(m);
	z=v_get(n);
	v=v_get(m);
	l=v_get(n);
	u=v_get(n);
	
	smrand(time(NULL));

	for(i=0;i<n;i++)
	{
		l->ve[i] = INT_MIN;
		u->ve[i] = INT_MAX;
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

	m_output(A);
	v_output(y);
	start = time(NULL);
	VEC *P = permutationreduction(A,y,l,u);
	stop=time(NULL);
	v_output(P);

	diff=difftime(stop,start);
	printf("Computation time is %e\n",diff);

	//V_FREE(P);
	M_FREE(A);
	V_FREE(z);
	V_FREE(v);
	V_FREE(y);
	
	return 0;
}
