/*
 * reduction.h
 *
 *  Created on: 2011-4-23
 *      Author: linwanru
 */

#ifndef REDUCTION_H_
#define REDUCTION_H_

void igt(MAT *A, MAT *Z, int i, int j);
MAT * reduction(MAT *R, MAT *Z, double delta, int n);
MAT * oldreduction(MAT *R, MAT *Z, double delta, int n);
#endif /* REDUCTION_H_ */
