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
	//In my matlab code, A is nxm, the library gives A mxn, so I'm switching it here
	int n = A->m;
	int m = A->n;
	int i,j,L;
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
	tempA = m_get(n,m);
	MAT *tempG = m_get(m,n);
	_m_copy(A,tempA,0,0);
	svd(tempA,U,V,diag);
	//Compute V*D*U'
	for (i=0;i<m;i++)
	{
		set_col(tempG,i,sv_mlt(1/diag->ve[i],get_row(V,i,VNULL),VNULL));
	}
	m_mlt(tempG,U,G);
	G = m_transp(G,MNULL);

	//Start the SW Algorithm
	VEC *tempVec = v_get(n);
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
			tempDbl = in_prod(y,get_col(G,i,VNULL));
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
			colNorm = v_norm2(get_col(G,i,VNULL));
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
		P->ve[L] = Index[best];
		tempInt = Index[L];
		Index[L] = Index[best];
		Index[best] = tempInt;
		sv_mlt(a,get_col(A,besti,VNULL),tempVec);
		y = v_sub(y,tempVec,VNULL);
		tempVec = sv_mlt(1/bestColNorm,get_col(G,besti,VNULL),VNULL);
		for (j=0;j<L;j++)
		{
			i = Index[j];
			tempDbl = in_prod(get_col(G,besti,VNULL),get_col(G,i,VNULL));
			set_col(G,i,v_sub(get_col(G,i,VNULL),sv_mlt(tempDbl,tempVec,VNULL),VNULL));
		}
	}
	
	V_FREE(tempVec);
	M_FREE(tempA);
	M_FREE(tempG);
	V_FREE(diag);
	M_FREE(U);
	M_FREE(V);	

	return P;
}
