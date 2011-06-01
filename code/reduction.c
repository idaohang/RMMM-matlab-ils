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
   set_col(Z,j,v_mltadd(get_col(Z,j,VNULL),get_col(Z,i,VNULL),-kesai,VNULL));
	}
}

MAT * reduction(MAT *R_y, MAT *Z, double delta, int n)
{

	int k=1;
	int i,j;
	Real c=0, s=0;

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
	    	   set_col(Z,k,v_mltadd(get_col(Z,k,VNULL),get_col(Z,k1,VNULL),-kesai,VNULL));

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
	for(j=0;j<n;j++)
	{
		if (R_y->me[j][j]<0)
			set_row(R_y,j,sv_mlt(-1, get_row(R_y,j,VNULL),VNULL));

		for (i=j+1;i<n;i++)
		{
			if(R_y->me[i][j]<100*MACHEPS )
			R_y->me[i][j]=0;
	    }
	}

	return(R_y);
	}
MAT * oldreduction(MAT *R, MAT *Z, double delta, int n)
{
	int k=1;
	int j;
	Real c=0, s=0;
	m_resize(R,n,n+1);

	while(k<n)
	{
		int k1=k-1;
		int i;
        int kesai=round(R->me[k1][k]/R->me[k1][k1]);
        igt(R,Z,k-1,k);
        /*Test if permutation is needed.*/
	    if (delta*R->me[k1][k1]*R->me[k1][k1] < R->me[k1][k]-kesai*R->me[k1][k1]+R->me[k][k]*R->me[k][k])
	    {
	    	if(kesai!=0){
	    		for(i=k-2;i>=0;i--)
	    		{
	    			igt(R,Z,i,k);
	    			/* update Z*/
	    			set_col(Z,k,v_mltadd(get_col(Z,k,VNULL),get_col(Z,k1,VNULL),-kesai,VNULL));
	    		}
	    		k=k+1;
	    	}
	    	else

	    	{
	    		pivot(R,k-1,k,1);
	    		pivot(Z,k-1,k,1);
	    		givens(R->me[k-1][k-1],R->me[k][k-1],&c,&s);
	    		rot_rows(R,k-1,k,c,s,R);
	    		 if (k>1){
	    			 k--;
	    			     }
	    	}

	}
/*   Interchange the rows of R to make all the diagnoal elements to be positive;*/
	for(j=0;j<n;j++)
	{
		if (R->me[j][j]<0)
			set_row(R,j,sv_mlt(-1, get_row(R,j,VNULL),VNULL));

		for (i=j+1;i<n;i++)
		{
			if(R->me[i][j]<100*MACHEPS )
			R->me[i][j]=0;
	    }
	  }
	}
	return(R);
	}
