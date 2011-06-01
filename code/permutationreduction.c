/*
 * permutationreduction.c
 * An implementation of the SW algorithm ("A New Ordering for Efficient Sphere Decoding"). There is a more efficient and stable
 * algorithm that gives the same result described in the paper by myself and Chang, 
 * but the code becomes more complex.
 */ 
#include <stdio.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "math.h"
#include "reduction.h"

/*
* Apply the SW Algorithm to the matrix A using the input vector y
* and boundary information defined in vectors l and u (lower and upper
* bound respectively)
*/
VEC * permutationreduction(MAT *A, VEC *y, VEC *l, VEC *u)
{
	int n = A->m;
	int m = A->n;
	int i,j,L,sig,a,best,count;
	double bestColNorm;	
	VEC *P = v_get(m);
	int Index[m];
	for (i =0;i<m;i++)
	{
		P->ve[i] = i;
		Index[i] = i;
	}
	
	//Compute G, the moore-penrose generalized inverse of A
	MAT *G = m_get(m,n);
	m_zero(G);
	MAT *tempA;
	VEC *diag;
	MAT *U,*V;
	U = m_get(n,n);
	V = m_get(m,m);
	diag = v_get(n);
	temp = m_get(n,m);
	MAT *tempG = m_get(m,n);
	_m_copy(A,tempA,0,0);
	svd(tempA,U,V,diag);
	//Compute V*D*U'
	for (i=0;i<m;i++)
	{
		set_col(tempG,i,sv_mlt(1/diag->ve[i],get_row(V,i,VNULL),VNULL));
	}
	m_mlt(tempG,U,G);

	for (L=m;L>=0;L--)
	{
		sig = -1;
		a = -1;
		best = -1;
		bestColNorm = 0;
		count = 0;
		for (j = 1;j<L;j++)
		{
			i = Index[j];
			count=count+1;
			
		}
	}
	
	
	V_FREE(Index);
	M_FREE(tempA);
	M_FREE(tempG);
	V_FREE(diag);
	M_FREE(U);
	M_FREE(V);	

	return P;
}
