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
	int i;	
	VEC *P = v_get(m);
	VEC *Index = v_get(m);
	for (i =0;i<m;i++)
	{
		P->ve[i] = i;
		Index->ve[i] = i;
	}
	
	//Compute G, the moore-penrose generalized inverse of A
	MAT *G = m_get(n,m);
	MAT *temp;
	VEC *diag;
	MAT *U,*V;
	U = m_get(n,n);
	V = m_get(m,m);
	diag = v_get(n);
	temp = m_get(m,n);
	_m_copy(A,temp,0,0);
	svd(temp,U,V,diag);
	//Compute V*D*U'
	VEC *curCol = v_get(n);
	for (i=0;i<m;i++)
	{
		set_col(V,i,sv_mlt(diag->ve[i],get_col(V,i,VNULL),VNULL));
	}
	m_mlt(V,U,G);
	printf("A:\n");
	m_output(A);
	printf("G:\n");
	m_output(G);
	
	
	V_FREE(Index);
	M_FREE(temp);
	V_FREE(diag);
	M_FREE(U);
	M_FREE(V);	

	return P;
}
