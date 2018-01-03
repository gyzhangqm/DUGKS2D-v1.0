#ifndef _ARK_D2Q9_H
#define _ARK_D2Q9_H

#include "ZeroReference.h"

int const

DV_Qu = 1,

DV_Qv = 9;

double const

MaxU = sqrt(1.5/Lambda0),

Eta = 0,

MaSpan = sqrt(3.0*RT);

#endif
// const double RT = 1.0/3.0,s_RT = sqrt(3.0*RT);
// //
// const typexi xi[Q]={{0,0},{s_RT,0},{s_RT,s_RT},{0,s_RT},
// 					{-s_RT,s_RT},{-s_RT,0},{-s_RT,-s_RT},{0,-s_RT},{s_RT,-s_RT}};
// //
// const double omega[Q]={4.0/9.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0};					
// #endif
