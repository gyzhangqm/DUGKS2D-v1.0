#ifndef _ZERO_REFERENCE_H
#define _ZERO_REFERENCE_H

#include <cmath>
#include "ZeroFlip.h"

double const

infinitesimal = 1.0E-14,

PI = 3.141592653589793;

//----------------------Mesh file---------------
const int NL = 16;//name of mesh file

const double

ChLength = 1,

MinL = 6.037849742228585e-02,//3.079505855617000e-02,//1.555181192035053e-02,
//7.776141069016656e-03,//1.0/NL,
X_Beg = 0, 

X_End = 1, 

Y_Beg = 0.0, 

Y_End = 1,

Lx = X_End - X_Beg,

Ly = Y_End - Y_Beg,

Kh3 = 1.0E-3*MinL*MinL*MinL; //Limiter

//------------------------------Atomic species-------------------------------------
const double 

Omega0 = 0.5, //hard sphere(HS) = 0.5, variable hard sphere(VHS) = 0.68

Pr = 2.0/3.0, //Prandtl Number

nK = 0,//internal degree 0 = single  2 = double

agEq = (nK + 3.0 - 2.0)/2.0,//used in the Equilibrium of g(x,xi,t)

Cv = (nK + 3.0)/2.0,

Gamma = (nK + 5.0)/(nK + 3.0);

const double 

T0 = 1.0,

R0 = 1.0,

RT = R0*T0,

Lambda0 = 1/(2.0*R0*T0),

Rho0 = 1.0,

U0 = 0.1,//W_i*TC_r,

Ma = U0/sqrt(R0*T0),

V0 = 0.0,

Re = 100,

Mu0 = Rho0*U0*ChLength/Re,

Nu0 = Mu0/Rho0,

Kn = 16.0*Mu0/(5.0*Rho0*R0*T0)*sqrt(1.0/(4.0*Lambda0*PI));

#endif
