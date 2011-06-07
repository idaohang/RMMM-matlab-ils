/*
 * qrmcp.c
 *
 */

#include <stdio.h>
#include <math.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include <time.h>
/*Find the index of minimum value among the k to last entities of a vector*/
u_int findmin(VEC* a, u_int k) {
	u_int mincol;
	u_int t;
	double tmp;
	tmp = a->ve[k];
	mincol = k;
	for (t = k + 1; t < a->dim; t++) {
		if (a->ve[t] < tmp) {
			tmp = a->ve[t];
			mincol = t;
		}

	}
	return mincol;

}
/* If type ==1,interchange the ith and jth column of Matrix A;*/
/*If type ==0,interchange the ith and jth row of Matrix A;*/
void pivot(MAT * A, int i, int j, int type) {
	VEC *tmp;
	if (type == 1) {
		tmp = get_col(A, j, VNULL);
		set_col(A,j,get_col(A,i,VNULL));
		set_col(A,i,tmp);
	} else if (type == 0) {
		tmp = get_row(A, j, VNULL);
		set_row(A,j,get_row(A,i,VNULL));
		set_row(A,i,tmp);
	}
	V_FREE(tmp);
}
/*Apply Householder factorization with minimum column pivoting on
matrix A,and y, Z is the permutation matrix. Return A=[R,y]=[Q^T*A,Q^T*y];
*/

MAT* qrmcp(MAT *A_y, MAT *Z,int n) {

	m_zero(Z);
	u_int i = 0;
	u_int k;
	u_int mincol;
	Real tao = 0;
	VEC *a_norm, *u;
	int* p = malloc(sizeof(int) * n);
	double tmp;

	u = v_get(A_y->m);
	a_norm = v_get(n);

	for (i = 0; i < n; i++)
		p[i] = i;

	/* compute the two norm of each column of A;*/
	for (i = 0; i < n; i++) {
		a_norm->ve[i] = v_norm2(get_col(A_y,i,VNULL))*v_norm2(get_col(A_y,i,VNULL));

	}

	for (i = 0; i < A_y->n - 1; i++) {
		/* Find the column with the minimum two norm and do pivoting on A and Z;*/
		mincol = findmin(a_norm, i);

		if (mincol > i) {
			pivot(A_y, i, mincol, 1);
			tmp = a_norm->ve[mincol];
			a_norm->ve[mincol] = a_norm->ve[i];
			a_norm->ve[i] = tmp;
			int pi = p[i];
			p[i] = p[mincol];
			p[mincol] = pi;

		}

		/*Compute the Householder vector u for the ith column of A;*/
		hhvec(get_col(A_y, i, VNULL), i, &tao, u, &A_y->me[i][i]);
		/*Apply Householder on A_y;*/
		hhtrcols(A_y, i, i + 1, u, tao);
		/* Zero the entities below diagonal;*/
		for (k = i + 1; k < A_y->m; k++) {
			A_y->me[k][i] = 0;
		}

		/*Compute the new norm*/
		for (k = i + 1; k < A_y->n - 1; k++) {
			a_norm->ve[k] =a_norm->ve[k] * a_norm->ve[k] - A_y->me[i][k]* A_y->me[i][k];

		}

	}
	/*Form the permutation matrix*/
	for (i = 0; i < n; i++) {
		Z->me[p[i]][i] = 1;
	}

	free(p);
	V_FREE(u);
	V_FREE(a_norm);
	return (A_y);
}

