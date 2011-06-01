/*
 * ILS.c
 *
 */
#include <stdio.h>
#include "meschach/matrix.h"
#include "qrmcp.h"
#include "reduction.h"
#include "search.h"
#include "string.h"
#include <time.h>
VEC * ils(MAT *A, VEC *y)
{
	MAT *Z, *R;
	VEC *z_tilt,*x;
	time_t t,start,stop;
	long double diff;
	start= time(NULL);
	int m, n;
	m = A->m;
	n = A->n;
    z_tilt=v_get(n);
	Z = m_get(n, n);
	R=m_get(m,n+1);

	_m_copy(A,R,0,0);
	set_col(R,n,y);
	x =v_get(n);
	/*	 Apply QR factorization with minimum column pivoting on A and y to get R.
	 R=[R,y]=[Q^T*A,Q^T*y]*/
	qrmcp(R,Z,n);
	m_resize(R,n,n+1);
    reduction(R, Z, 1, n);
	search(R, n,z_tilt);
	mv_mlt(Z, z_tilt, x);
	M_FREE(Z);
	M_FREE(R);
    V_FREE(z_tilt);
    stop=time(NULL);
    diff=difftime(stop,start);
    printf("Computation time is %e\n",diff);
    return(x);
}
