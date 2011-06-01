/*
 * qrmcp.h
 *
 *  Created on: 2011-4-22
 *      Author: linwanru
 */

#ifndef QRMCP_H_
#define QRMCP_H_
u_int findmin(VEC* a,int k);
void pivot(MAT * A,int i, int j, int type);
MAT* qrmcp(MAT *A_y, MAT *Z,int n);

#endif /* QRMCP_H_ */

