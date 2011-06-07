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
PERM * permutationreduction(MAT *A, VEC *y, VEC *l, VEC *u)
{
	/*In my matlab code, A is nxm, the library gives A mxn, so I'm switching it here*/
	int n = A->m;
	int m = A->n;
	int i,j,L;
	PERM *P = px_get(m);
	PERM *Pinv = px_get(m);
	px_ident(P);
	int Index[m];
	for (i =0;i<m;i++)
	{
		Index[i] = i;
	}
	/*Compute G, the moore-penrose generalized inverse of A*/
	MAT *G = m_get(m,n);
	MAT *tempA;
	VEC *diag;
	MAT *U,*V;
	U = m_get(n,n);
	V = m_get(m,m);
	diag = v_get(n);
	tempA = m_get(n,m);
	MAT *tempG = m_get(m,n);
	VEC *tempVecN = v_get(n);
	VEC *tempVecN2 = v_get(n);
	VEC *tempVecM = v_get(m);
	VEC *tempVecN3 = v_get(n);

	m_copy(A,tempA);
	svd(tempA,U,V,diag);
	/*Compute V*D*U'*/
	for (i=0;i<m;i++)
	{
		get_row(V,i,tempVecM);
		sv_mlt(1/diag->ve[i],tempVecM,tempVecM);
		set_col(tempG,i,tempVecM);
	}
	m_mlt(tempG,U,G);
	m_transp(G,tempG);

	/*Start the SW Algorithm*/
	int best,besti,count,tempInt,a,ap,bp;
	double tempDbl,colNorm,sigp,sig,bestColNorm;
	for (L=m-1;L>=0;L--)
	{
		sig = -1;
		a = -1;
		best = -1;
		bestColNorm = 0;
		count = -1;
		for (j = 0;j<=L;j++)
		{
			i = Index[j];
			count=count+1;
			get_col(G,i,tempVecM);
			tempDbl = in_prod(y,tempVecM);
			ap = max(min((int)round(tempDbl),(int)u->ve[i]),(int)l->ve[i]);
		
			if(ap == l->ve[i])
			{
				bp = ap+1;
			}
			else if(ap == u->ve[i])
			{
				bp = ap-1;
			}
			else
			{
				if((tempDbl - (double)ap) < 0)
				{
					tempInt = -1;
				}
				else
				{
					tempInt = 1;
				}
				bp = ap + tempInt;
			}
			get_col(tempG,i,tempVecM);
			colNorm = v_norm2(tempVecM);
			sigp = (1/colNorm)*(fabs(tempDbl-(double)bp));
			if(sigp >= sig)
			{
				a = ap;
				sig = sigp;
				best = j;
				besti=i;
				bestColNorm = colNorm*colNorm;
			}
		}
		px_transp(P,L,best);
		tempInt = Index[L];
		Index[L] = Index[best];
		Index[best] = tempInt;
		get_col(A,besti,tempVecN);
		sv_mlt(a,tempVecN,tempVecN);
		v_sub(y,tempVecN,tempVecN2);
		v_copy(y,tempVecN2);
		get_col(tempG,besti,tempVecN);
		sv_mlt(1/bestColNorm,tempVecN,tempVecN);
		for (j=0;j<L;j++)
		{
			i = Index[j];
			get_col(tempG,besti,tempVecN);
			get_col(tempG,i,tempVecN2);
			tempDbl = in_prod(tempVecN,tempVecN2);
			
			sv_mlt(tempDbl,tempVecN,tempVecN);
			get_col(tempG,i,tempVecN2);
			v_sub(tempVecN2,tempVecN,tempVecN3);
			set_col(tempG,i,tempVecN3);
		}
	}
	
	px_inv(P,Pinv);
	V_FREE(tempVecN);
	M_FREE(tempA);
	M_FREE(tempG);
	V_FREE(diag);
	M_FREE(U);
	M_FREE(V);
	M_FREE(G);	
	V_FREE(tempVecM);
	V_FREE(tempVecN2);
	V_FREE(tempVecN3);

	return Pinv;
}
