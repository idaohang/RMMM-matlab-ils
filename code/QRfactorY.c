#include	<stdio.h>
#include	<math.h>
#include        "meschach/matrix2.h"

/* QRfactor -- forms the QR factorisation of A -- factorisation stored in
   compact form as described above ( not quite standard format ) */
MAT	*QRfactorY(MAT *A, VEC *diag, VEC *y)
{
    unsigned int	k,limit;
    Real	beta;
    static	VEC	*hh=VNULL, *w=VNULL;
    
    if ( ! A || ! diag )
	error(E_NULL,"QRfactor");
    limit = min(A->m,A->n);
    if ( diag->dim < limit )
	error(E_SIZES,"QRfactor");
    
    hh = v_resize(hh,A->m);
    w  = v_resize(w, A->n);
    MEM_STAT_REG(hh,TYPE_VEC);
    MEM_STAT_REG(w, TYPE_VEC);
    
    for ( k=0; k<limit; k++ )
    {
	/* get H/holder vector for the k-th column */
	get_col(A,k,hh);
	/* hhvec(hh,k,&beta->ve[k],hh,&A->me[k][k]); */
	hhvec(hh,k,&beta,hh,&A->me[k][k]);
	diag->ve[k] = hh->ve[k];
	
	/* apply H/holder vector to remaining columns */
	/* hhtrcols(A,k,k+1,hh,beta->ve[k]); */
	hhtrcols(A,k,k+1,hh,beta);
	hhtrvec(hh,beta,k,y,y);
    }

#ifdef	THREADSAFE
    V_FREE(hh);	V_FREE(w);
#endif

    return (A);
}
