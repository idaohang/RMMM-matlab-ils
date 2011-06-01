/*
*randn.c
*Use the polar form of the Box-Muller transform to generate random numbers
*from a 0-mean standard deviation 1 normal distribution
*/ 
#include <stdio.h>
#include <time.h>
#include "math.h"
#include <stdlib.h>

double randn()
{ 
	float x1, x2, w, y1;
	do
	{
		x1 = 2.0 * (double)rand()/RAND_MAX - 1.0;
		x2 = 2.0 * (double)rand()/RAND_MAX - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	return y1;
}
