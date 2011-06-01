/*
 * search.c
 */
# include <stdio.h>
# include <math.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include <time.h>
int sgn(double a) {
	if (a > 0)
		return 1;
	else
		return -1;
}

VEC * search(MAT * R_y,  int n, VEC *z_hat) {

	MAT * S;
	VEC  *c, *prsd;
	IVEC *d,*z;
	double beta = __builtin_inf(), gama;
	int k = n - 1, i = 0;

	S = m_get(n, n + 1);
	z = iv_get(n);
	c = v_get(n);
	d = iv_get(n);
	prsd = v_get(n);
	if(z_hat==VNULL)
		z_hat = v_get(n);

	m_zero(S);
	iv_zero(z);
	v_zero(c);
	iv_zero(d);
	v_zero(prsd);

	double newsprsd;

	c->ve[k] = R_y->me[k][n] / R_y->me[k][k];
	z->ive[k] = round(c->ve[k]);
	gama = R_y->me[k][k] * (c->ve[k] - z->ive[k]);
	d->ive[k] = sgn(c->ve[k] - z->ive[k]);

	while (TRUE) {

		newsprsd = prsd->ve[k] + gama * gama;
		if (newsprsd < beta) {
			if (k != 0) {
				for (i = 0; i <= k; i++) {
					S->me[i][k] = R_y->me[i][k] * z->ive[k] + S->me[i][k + 1];
				}
				k--;
				prsd->ve[k] = newsprsd;
				c->ve[k] = (R_y->me[k][n] - S->me[k][k + 1]) / R_y->me[k][k];
				z->ive[k] = round(c->ve[k]);
				gama = R_y->me[k][k] * (c->ve[k] - z->ive[k]);
				d->ive[k] = sgn(c->ve[k] - z->ive[k]);

			} else {
				int index;
				for(index=0;index<z->dim;index++)
					z_hat->ve[index]=z->ive[index];
				beta = newsprsd;
				k = 1;
				z->ive[k] = z->ive[k] + d->ive[k];
				gama = R_y->me[k][k] * (c->ve[k] - z->ive[k]);
				d->ive[k] = -d->ive[k] - sgn(d->ive[k]);
			}
		} else {
			if (k == n - 1)
				break;
			else {
				k++;
				z->ive[k] = z->ive[k] + d->ive[k];
				gama = R_y->me[k][k] * (c->ve[k] - z->ive[k]);
				d->ive[k] = -d->ive[k] - sgn(d->ive[k]);

			}
		}
	}

	M_FREE(S);
	IV_FREE(z);
	V_FREE(c);
	IV_FREE(d);
	V_FREE(prsd);

	return z_hat;
}
