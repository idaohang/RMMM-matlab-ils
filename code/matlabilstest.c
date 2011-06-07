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
#include <limits.h>
#include "qrmcp.h"
#include "reduction.h"
#include "search.h"
#include "mex.h"
#include "mat.h"

void VEC2mxArray(VEC *in, const mxArray *out)
  {
    int rows,i;
    double* temp = (double *) mxGetPr(out);
    if (temp==0) mexErrMsgTxt("vec2mxArray: Pointer to data is NULL");
    rows = in->dim;
    if (rows==0) mexErrMsgTxt("vec2mxArray: Data has zero rows");
	for(i=0; i<rows; i++) {
	  *temp++= in->ve[i];
	}
  }

void MAT2mxArray(MAT *in, const mxArray *out)
  {
    int rows, cols, i, j;

    double* temp = (double *) mxGetPr(out);
    if (temp==0) mexErrMsgTxt("mat2mxArray: Pointer to data is NULL");

    rows = in->m;
    cols = in->n;
    if (rows==0) mexErrMsgTxt("mat2mxArray: Data has zero rows");
    if (cols==0) mexErrMsgTxt("mat2mxArray: Data has zero columns");

    for(j=0; j<cols; j++) {
	for(i=0; i<rows; i++) {
	  *temp++= in->me[i][j];
	}
    }

  }
MAT *mxArray2MAT(const mxArray *in) {
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

VEC * mxArray2VEC(const mxArray *in)
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

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	MAT *A;
	VEC *y,*low,*up;
	clock_t start,stop;
	double ratio,diff;
	ratio = 1.0/CLOCKS_PER_SEC;
	
	int m,n,i,j;
	mwSize rows,cols;

	/*Convert matlab data to meschach*/
	A = mxArray2MAT(prhs[0]);
	y = mxArray2VEC(prhs[1]);

	m= A->m;
	n=A->n;

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
	VEC *z_ils = v_get(n);
	VEC *diag = v_get(n);
	VEC *yp = v_get(m);
	MAT * Rp = m_get(n,n);
	
	_m_copy(A,R,0,0);
	set_col(R,n,y);
	
	start = clock();
	/*	Apply QR factorization with minimum column pivoting on A and y to get R.
	 R=[R,y]=[Q^T*A,Q^T*y]*/
	qrmcp(R,Z,n);
	m_resize(R,n,n+1);
    reduction(R, Z, 1, n);
	get_col(R,n,yp);
	m_resize(R,n,n);

	/*Find the permutations, permute R, restore to upper-triangular & apply hh vecs to y*/
	PERM *P = permutationreduction(R,yp,low,up);
	px_cols(P,R,Rp);
	QRfactorY(Rp,diag,yp);
	for (i=0;i<n;i++)
	{
		for(j=i+1;j<n;j++)
		{
			Rp->me[j][i] = 0;
		}
	}

	/*Append y to R again*/
	m_resize(Rp,n,n+1);
	set_col(Rp,n,yp);
	/*Search*/
	search(Rp,n,z_tilt);
	/*Transform solution*/
	px_vec(P,z_tilt,x);
	mv_mlt(Z,x,z_ils);
	stop=clock();
	
	diff=ratio*(long)stop - ratio*(long)start;
	
	plhs[0] = mxCreateDoubleScalar(diff);
	cols=1;
	rows=n;
	plhs[1] = mxCreateDoubleMatrix(rows,cols,mxREAL);
	VEC2mxArray(z_ils,plhs[1]);
	
	PX_FREE(P);
	M_FREE(A);
	M_FREE(Rp);
	V_FREE(low);
	V_FREE(up);
	V_FREE(diag);
	V_FREE(z_ils);
	V_FREE(y);
	M_FREE(R);
	M_FREE(Z);
	V_FREE(z_tilt);
	V_FREE(x);
	V_FREE(yp);
}
