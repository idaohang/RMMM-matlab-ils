/*
 * test.c
 *
 */
#include <stdio.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "math.h"
#include "ils.h"
#include "mils.h"
#include <time.h>




void main()
{
        MAT  *A,*B;
		VEC *z, *x,*y;
		time_t t,start,stop;
		double diff;
		int i,m,k,n;
        start= time(NULL);
        m=20;
        k=10;n=10;

		A=m_get(m,k);
		B=m_get(m,n);
		y=v_get(m);
		x=v_get(k);
		z=v_get(n);

		smrand((unsigned)time(&t));
        m_rand(A);
        m_rand(B);
        v_rand(x);
        v_rand(z);
        for(i=0;i<n;i++)
        {
        	z->ve[i]=round(10*z->ve[i]);
        }
        v_add(mv_mlt(A,x,VNULL), mv_mlt(B,z,VNULL), y);

        printf("For y=A*x+B*Z,the data used for test:\n");
        printf("A is:\n");
        m_output(A);
        printf("B is:\n");
        m_output(B);
        printf("y is:\n");
        v_output(y);
        printf("x used for test is:\n");
        v_output(x);
        printf("z used for test is:\n");
        v_output(z);

	    printf("Miles solution is\n");

	   	v_output(mils(A,B,y));

	    M_FREE(A);
	    M_FREE(B);
	    V_FREE(z);
	    V_FREE(x);
	    V_FREE(y);
        stop=time(NULL);
        diff=difftime(stop,start);
        printf("Computation time is %e\n",diff);


}
