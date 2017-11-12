#ifndef _ZERO_REFERENCE_H
#define _ZERO_REFERENCE_H

#include <cmath>
#include "ZeroFlip.h"

double const

infinitesimal = 1.0E-14,

PI = 3.14159265358979;

//----------------------Mesh file---------------
const int NL = 100;//name of mesh file

const double

ChLength = 1,

MinL = ChLength/NL,//4.999999999999893E-1,//1.0/NL,

X_Beg = -0.5, 

X_End = 0.5, 

Y_Beg = 0.0, 

Y_End = 0.02,

Lx = X_End - X_Beg,

Ly = Y_End - Y_Beg,

Kh3 = 1.0E-3*MinL*MinL*MinL; //Limiter

//------------------------------Atomic species-------------------------------------
const double 

Omega0 = 0.5, //hard sphere(HS) = 0.5, variable hard sphere(VHS) = 0.68

Pr = 2.0/3.0, //Prandtl Number

nK = 2,//internal degree 0 = single  2 = double

Cv = (nK + 3.0)/2.0,

Gamma = (nK + 5.0)/(nK + 3.0);

const double 

T0 = 1.0,

R0 = 1.0,

Lambda0 = 1/(2.0*R0*T0),

Rho0 = 1.0,

Ma = 0.0,

U0 = Ma*sqrt(Gamma*R0*T0),//W_i*TC_r,

V0 = 0.0,

Mu0 = 1E-4,

Nu0 = Mu0/Rho0,

Re = Rho0*U0*ChLength/Mu0,

Kn = 16.0*Mu0/(5.0*Rho0*R0*T0)*sqrt(1.0/(4.0*Lambda0*PI));

#endif