#ifndef _ZERO_INFORMATION_H_
#define _ZERO_INFORMATION_H_

#include "D2GHn.h"

//------------------------------Normalized Parameters----------------------------
const double 

CFL = 0.1,

dt = CFL > 0.0 ? CFL*MinL/(MaxU): 1.0E-4,

h = 0.5*dt;

//---------------------Rankine-Hugoniot-----------------
/*const double

Ma_Outlet = sqrt((Ma*Ma*(Gamma-1.0)+2.0)/(2.0*Gamma*Ma*Ma-(Gamma-1))),

Rho_Outlet = Rho0*(((Gamma+1)*Ma*Ma)/((Gamma-1)*Ma*Ma+2)),

T_Outlet = T0*((1+0.5*(Gamma-1)*Ma*Ma)*(2.0*Gamma/(Gamma-1)*Ma*Ma-1)/
			(Ma*Ma*(2.0*Gamma/(Gamma-1)+0.5*(Gamma-1)))),

U_Outlet = Ma_Outlet*sqrt(Gamma*R0*T_Outlet),

V_Outlet = 0.0;*/

//---------------------2D-Riemann---------------------------------

const double

Rho1 = 0.5313, U1 = 0.0000,V1 = 0.0000,p1 = 0.4,T1 = p1/(Rho1*R0),

Rho2 = 1.0000, U2 = 0.7276,V2 = 0.0000,p2 = 1.0,T2 = p2/(Rho2*R0),

Rho3 = 0.8000, U3 = 0.0000,V3 = 0.0000,p3 = 1.0,T3 = p3/(Rho3*R0),

Rho4 = 1.0000, U4 = 0.0000,V4 = 0.7276,p4 = 1.0,T4 = p4/(Rho4*R0);

//-------------------------Shock Tube---------------------------

const double 

Rho_Outlet = 0.125,

U_Outlet = 0.0,

V_Outlet = 0.0,

T_Outlet = 0.8;

/*const int

I_TimeEnd = log(2.0)/(8.0*PI*PI*nu*dt);//TimeEnd/dt;

const double

TimeEnd = I_TimeEnd*dt;*/
//------------------------------Taylor-Couette--------------------------------
const double 

TC_R = 1.5,

TC_r = 0.5,

W_o = 0.0,

W_i = 0.2;

//-----------------------------Output-------------------
const int

VelocityZone = 7,//7 == TC

End_Step = 5000000,//log(2.0)/(8.0*PI*PI*Nu0*dt),

ZeroDebugControl = 1000, //

ConvergenceControl = 100, //

ResidualControl = 100, //print to screen

writeFileControl = 10000;

double const

RESIDUAL = 1.0E-8;

#endif
