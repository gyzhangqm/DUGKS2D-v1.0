#ifndef _D2Q9_ARK_H
#define _D2Q9_ARK_H
const int Q = 9;
#include <cmath>
struct typexi
{
	typexi(double a,double b):u(a),v(b){}
	double  u;
	double  v;
};
const double RT = 1.0/3.0,s_RT = sqrt(3.0*RT);
//
const typexi xi[Q]={{0,0},{s_RT,0},{s_RT,s_RT},{0,s_RT},
					{-s_RT,s_RT},{-s_RT,0},{-s_RT,-s_RT},{0,-s_RT},{s_RT,-s_RT}};
//
const double omega[Q]={4.0/9.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0};					
#endif
