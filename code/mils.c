/*
 * mils.c
 *
 */
#include <stdio.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "ils.h"
#include "search.h"

VEC * mils(MAT *A, MAT *B, VEC *y) {
	Real tao = 0;
	double diag = 0.0;
    int m,n,k,i,j;
    VEC *y_km,*z_hat,*x,*u;
    MAT *B_km;
    m=A->m;
    k=A->n;
    n=B->n;

    x=v_get(k);
    y_km=v_get(m-k);
    z_hat=v_get(n);
    B_km=m_get(m-k,n);
    u=v_get(m);
     /*Apply houshold transformation on A to get the QR factorization of A. Form Q^T*y and Q^T*B*/
    for(j=0;j<k;j++){

    	hhvec(get_col(A,j,VNULL),j,&tao,u,&A->me[j][j]);
    	hhtrcols(A,j,j+1,u,tao);
    	for(i=j+1;i<m;i++)
    	{
    		A->me[i][j]=0;
    	}

    	y=hhtrvec(u,tao,j,y,VNULL);
    	hhtrcols(B,j,0,u,tao);
    }

    v_move(y,k,m-k,y_km,0);
    v_resize(y,k);

    m_move (B,k,0,m-k,n,B_km,0,0);
    m_resize(B,k,n);
    z_hat=ils(B_km,y_km);
    m_resize(A,k,k);
    Usolve(A, v_sub(y, mv_mlt(B, z_hat, VNULL), VNULL),x, diag);
    v_resize(x,n+k);
    for(i=k;i<n+k;i++)
    {
    	x->ve[i]=z_hat->ve[i-k];
    }
    V_FREE(y_km);
    M_FREE(B_km);
    V_FREE(z_hat);
    return x;
}


