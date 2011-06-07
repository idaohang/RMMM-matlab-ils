/*
 * reduction.c
 *
 */

# include <stdio.h>
# include <math.h>
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "qrmcp.h"
#include <time.h>
Real  macheps = MACHEPS;
VEC *tmp1,*tmp2,*tmp3;
void igt(MAT *A, MAT *Z, int i, int j)
{
	int kesai,k;
	kesai=round(A->me[i][j]/A->me[i][i]);
	if(kesai!=0){
/*   Subtract the jth column by kesai times ith colum;*/
	for(k=0;k<=i;k++)
	{
		A->me[k][j]=A->me[k][j]- kesai*A->me[k][i];
	}
 /* update Z*/
	get_col(Z,j,tmp1);
	get_col(Z,i,tmp2);
	v_mltadd(tmp1,tmp2,-kesai,tmp3);
   set_col(Z,j,tmp3);
	}
}

MAT * reduction(MAT *R_y, MAT *Z, double delta, int n)
{

	int k=1;
	int i,j;
	Real c=0, s=0;
	tmp1 = v_get(n);
	tmp2 = v_get(n);
	tmp3 = v_get(n);

	while(k<n)
	{
		int k1=k-1;
		int i;
        int kesai=round(R_y->me[k1][k]/R_y->me[k1][k1]);
        double t=R_y->me[k1][k]-kesai*R_y->me[k1][k1];
        /*Test if permutation is needed.*/
	    if (delta*R_y->me[k1][k1]*R_y->me[k1][k1] > t*t+R_y->me[k][k]*R_y->me[k][k])
	    {
	    	if(kesai!=0){
	    		for(i=0;i<k1;i++)
	    		{
	    			R_y->me[i][k] -= kesai*R_y->me[i][k1];
	    		}
	    		R_y->me[k1][k]=t;
	    	  /* update Z*/
	    	  get_col(Z,k,tmp1);
	    	  get_col(Z,k1,tmp2);
	    	  v_mltadd(tmp1,tmp2,-kesai,tmp3);
	    	   set_col(Z,k,tmp3);

	    	   if(kesai>=2 || kesai<=-2)
	    	   {
	    		for (i=k-2;i>=0;i--)
	    		{
	    			igt(R_y,Z,i,k);

               }
	    	   }
	    	}
			 pivot(R_y,k-1,k,1);
		     givens(R_y->me[k-1][k-1],R_y->me[k][k-1],&c,&s);
		     rot_rows(R_y,k-1,k,c,s,R_y);

		     pivot(Z,k-1,k,1);
	         if (k>1)
	        	 k--;
		     }
	   else
	   {
		   k++;
	   }
	}
	
	VEC * tmp4 = v_get(R_y->n);
	VEC * tmp5 = v_get(R_y->n);
	for(j=0;j<n;j++)
	{
		if (R_y->me[j][j]<0){
			get_row(R_y,j,tmp4);
			sv_mlt(-1, tmp4,tmp5);
			set_row(R_y,j,tmp5);
		}

		for (i=j+1;i<n;i++)
		{
			if(R_y->me[i][j]<100*MACHEPS )
			R_y->me[i][j]=0;
	    }
	}
	
	V_FREE(tmp1);
	V_FREE(tmp2);
	V_FREE(tmp3);
	V_FREE(tmp4);
	V_FREE(tmp5);
	return(R_y);
	}
