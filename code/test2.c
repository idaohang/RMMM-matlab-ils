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
#include "permutationreduction.h"

int main()
{
      MAT *A;
		VEC *z,*y,*v,*l,*u;
		time_t t,start,stop;
		double diff;
		int m,n,i;
		double sigma = 0.5;
		m=5;
		n=3;

		srand(time(NULL));
		start = time(NULL);
		A=m_get(m,n);
		y=v_get(m);
		z=v_get(n);
		v=v_get(m);
		l=v_get(n);
		u=v_get(n);
		
		smrand((unsigned)time(&t));

		m_rand(A);
		v_rand(z);
		for(i=0;i<n;i++)
		{
			z->ve[i]=round(10*z->ve[i]);
			l->ve[i] = -10;
			u->ve[i] = 10;
		}
	
		for(i=0;i<m;i++)
		{
			v->ve[i]=randn()*sigma;
		}


		v_add(mv_mlt(A,z,VNULL), v, y);

/*		printf("A is:\n");
		m_output(A);
		printf("y is:\n");
		v_output(y);
		printf("z used for test is:\n");
		v_output(z);
		printf("noise is:\n");
		v_output(v);
*/
		VEC *P = permutationreduction(A,y,l,u);
		v_output(P);

		V_FREE(P);
		M_FREE(A);
		V_FREE(z);
		V_FREE(v);
		V_FREE(y);

		stop=time(NULL);
		diff=difftime(stop,start);
		printf("Computation time is %e\n",diff);

		return 0;
}
