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
#include "meschach/matlab.h"
#include "permutationreduction.h"
#include <limits.h>
#include "qrmcp.h"
#include "reduction.h"
#include "search.h"
#include "mex.h"
#include "mat.h"

MAT *mxArray2MAT(const mxArray *in)
{
	MAT *out;
	int i, j, rows, cols;
	double *temp = (double *) mxGetPr(in);

	if (temp == 0)
	mexErrMsgTxt("mxArray2mat: Pointer to data is NULL");
	rows = (int) mxGetM(in);

	if (rows == 0)
	mexErrMsgTxt("mxArray2mat: Data has zero rows");

	cols = (int) mxGetN(in);
	if (cols == 0)
	mexErrMsgTxt("mxArray2mat: Data has zero columns");
	out=m_get(rows,cols);
	for (i = 0; i < cols; i++) {
		for (j = 0; j < rows; j++) {
			out->me[j][i] = (*temp++);
		}
	}
	return out;
}

VEC *mxArray2VEC(const mxArray *in)
{
 VEC *out;
 int i, m;
 double *temp = (double*) mxGetPr(in);
 if (temp==0) mexErrMsgTxt("mxArray2vec: Pointer to data is NULL");
 m = (int)mxGetNumberOfElements(in);

 if (m==0) mexErrMsgTxt("mxArray2vec: Size of data is zero");

 out=v_get(m);

 for (i=0; i<m; i++) {
 	out->ve[i] = (*(temp++));
 }
 return out;
}

void mexfunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	MAT *A;
	VEC *z,*y,*v,*low,*up;
	time_t t,start,stop;
	double diff;
	int m,n,i,j;
	/*Convert matlab data to meschach*/
	A = mxArray2MAT(prhs[0]);
	y = mxArray2VEC(prhs[1]);

	m= A->m;
	n=A->n;

	srand(time(NULL));
	smrand(time(NULL));
	low=v_get(n);
	up=v_get(n);
	
	for(i=0;i<n;i++)
	{
		low->ve[i] = INT_MIN;
		up->ve[i] = INT_MAX;
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
	/*	Apply QR factorization with minimum column pivoting on A and y to get R.
	 R=[R,y]=[Q^T*A,Q^T*y]*/
	qrmcp(R,Z,n);
	m_resize(R,n,n+1);
   reduction(R, Z, 1, n);
 
   VEC *yp = v_get(m);
	yp = get_col(R,n,VNULL);
	m_resize(R,n,n);

	/*Find the permutations, permute R, restore to upper-triangular & apply hh vecs to y*/
	PERM *P = permutationreduction(R,yp,low,up);
	R = px_cols(P,R,MNULL);
	VEC *diag = v_get(n);
	QRfactorY(R,diag,yp);

	/*Append y to R again*/
	m_resize(R,n,n+1);
	set_col(R,n,yp);
	/*Search*/
	search(R,n,z_tilt);
	P = px_inv(P,PNULL);
	z_tilt = px_vec(P,z_tilt,VNULL);
	mv_mlt(Z,z_tilt,x);
	stop=time(NULL);
	diff=difftime(stop,start);

	plhs[0] = mxCreateDoubleScalar(diff);		
	
	M_FREE(A);
	V_FREE(z);
	V_FREE(v);
	V_FREE(y);
	M_FREE(R);
	M_FREE(Z);
	V_FREE(z_tilt);
	V_FREE(x);
	V_FREE(yp);
}
