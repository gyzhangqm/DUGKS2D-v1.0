#include <iostream>
#include "DUGKSDeclaration.h"
using std::cout;

double const degrees = 3.0 + nK;//total molecule degrees of freedom;

double const degPlus2 = degrees + 2.0;

namespace D2V16{
double const c = MaSpan,e = Eta;
double const c2 = c*c,e2 = e*e,c2e2 = c2*e2;
const double
	A = 1.0/(4.0*e2),
	B = (-21.0*c2 - e2)/(48.0*c2e2),
	D = (7.0*c2 + e2)/(32.0*c2e2),
	E = (-3.0*c2 - e2)/(96.0*c2e2),
	F = 1.0/(3.0*c2),
	G = (3.0*c2 - 5.0*e2)/(4.0*e2*(e2-3.0*c2)),
	H = 3.0*(c2 + e2)/(4.0*e2*(e2 - 3.0*c2)),
	I = (21.0*c2 - 11.0*e2)/(48.0*c2e2),
	J = (21.0*c2*c2 - 32.0*c2e2 + 11.0*e2*e2)/(32.0*c2e2*(e2-3.0*c2)),
	K = (21.0*c2*c2 - 36.0*c2e2 + 7.0*e2*e2)/(32.0*c2e2*(e2-3.0*c2)),
	L = (3.0*c2 - 5.0*e2)/(96.0*c2e2),
	N = 1.0/(2.0*c*e2),
	P = (3.0*c2 - e2)/(12.0*c2e2*c),
	Q = (c2 + e2)/(16.0*c2e2*c),
	R = (c2 - e2)/(16.0*c2e2*c),
	S = 1.0/(24.0*c2*c2),
	T = 1.0/(4.0*c2*(e2 - 3.0*c2)),
	U = (c2 - e2)/(32.0*c2*c2*(e2 - 3.0*c2)),
	V = (5.0*c2 - e2)/(32.0*c2*c2*(e2 - 3.0*c2));
//
const double iArkC [DV_Qv][DV_Qv] = 
{
	{0,0,0,A,0,G,H,N,0,-N,0,0,-N,0,T,-T},
	{0,0,0,A,0,H,G,0,N,0,-N,-N,0,0,-T,T},
	{0,0,0,A,0,G,H,-N,0,N,0,0,N,0,T,-T},
	{0,0,0,A,0,H,G,0,-N,0,N,N,0,0,-T,T},
	{2.0/3.0, 1.0/(3.0*c), 1.0/(3.0*c), B,F,I,I,-N/2,-N/2,P,P,N/2,N/2,-S,S/2,S/2},
	{2.0/3.0,-1.0/(3.0*c),1.0/(3.0*c),B,-F,I,I,N/2,-N/2,-P,P,N/2,-N/2,S,S/2,S/2},
	{2.0/3.0,-1.0/(3.0*c),-1.0/(3.0*c),B,F,I,I,N/2,N/2,-P,-P,-N/2,-N/2,-S,S/2,S/2},
	{2.0/3.0,1.0/(3.0*c), -1.0/(3.0*c),B,-F,I,I,-N/2,N/2,P,-P,-N/2,N/2,S,S/2,S/2},
	{-0.5,0,0,D,0,J,K,-N/8,0,Q,0,0,R,0,U,V},
	{-0.5,0,0,D,0,K,J,0,-N/8,0,Q,R,0,0,V,U},
	{-0.5,0,0,D,0,J,K,N/8,0,-Q,0,0,-R,0,U,V},
	{-0.5,0,0,D,0,K,J,0,N/8,0,-Q,-R,0,0,V,U},
	{1.0/12,-1.0/(24*c),-1.0/(24*c),E,-F/16,L,L,N/16,N/16,-P/8,-P/8,-R/2,-R/2,S/4,S/4,S/4},
	{1.0/12,1.0/(24*c),-1.0/(24*c),E,F/16,L,L,-N/16,N/16,P/8,-P/8,-R/2,R/2,-S/4,S/4,S/4},
	{1.0/12,1.0/(24*c),1.0/(24*c),E,-F/16,L,L,-N/16,-N/16,P/8,P/8,R/2,R/2,S/4,S/4,S/4},
	{1.0/12,-1.0/(24*c),1.0/(24*c),E,F/16,L,L,N/16,-N/16,-P/8,P/8,R/2,-R/2,-S/4,S/4,S/4}
};
const double xi_u[DV_Qv] = 
{
	c,0,-c,0,
	c,-c,-c,c,
	2*c,0,-2*c,0,
	2*c,-2*c,-2*c,2*c
};
const double xi_v[DV_Qv] = 
{
	0,c,0,-c,
	c,c,-c,-c,
	0,2*c,0,-2*c,
	2*c,2*c,-2*c,-2*c
};
const double xi2Eta2[DV_Qv] = 
{
	c2+e2,c2+e2,c2+e2,c2+e2,
	2*c2,2*c2,2*c2,2*c2,
	4*c2,4*c2,4*c2,4*c2,
	8*c2,8*c2,8*c2,8*c2
};
}

using D2V16::iArkC;

using D2V16::xi2Eta2;

extern double * const xi_u = new double[DV_Qv];

extern double * const xi_v = new double[DV_Qv];

//-----------------------------functions to be defined----------------------------

/*
void DiscreteVelocityAssign()

void setXiDotdS()

void Update_phi_Eq(Cell_2D &cell)

void Update_phi_Eqh(Face_2D &face)

void Update_MacroVar(Cell_2D& cell)

void Update_MacroVar_h(Face_2D& face)

void Update_phi_fBh(Face_2D &face)

void Update_phi_fBh_Limiter(Face_2D &face)*/



void DiscreteVelocityAssign()
{
	for(int j = 0;j < DV_Qv;++j)
	{
		xi_u[j] = D2V16::xi_u[j];
		xi_v[j] = D2V16::xi_v[j];
	}
}

void setXiDotdS()
{
	for(int n = 0;n < Faces;++n)
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		FaceArray[n].xi_n_dS[i][j] = 
		FaceArray[n].Area*(xi_u[j]*FaceArray[n].Vx + xi_v[j]*FaceArray[n].Vy);
	}
}
void Macro16(double *M16,double const &rho,double const &u,double const &v,
				double const &T,double const &p)
{
	double const u2 = u*u, v2 = v*v, RT = R0*T;
	M16[0]  =  rho;
	M16[1]  =  rho*u;
	M16[2]  =  rho*v;
	M16[3]  =  rho*(degrees*RT + u2 + v2);
	M16[4]  =  M16[1]*v;
	M16[5]  =  rho*u2 + p;
	M16[6]  =  rho*v2 + p;
	M16[7]  =  M16[1]*(degPlus2*RT + u2 + v2);
	M16[8]  =  M16[2]*(degPlus2*RT + u2 + v2);
	M16[9]  =  M16[1]*(3.0*RT + u2);
	M16[10] =  M16[2]*(3.0*RT + v2);
	M16[11] =  M16[2]*(RT + u2);
	M16[12] =  M16[1]*(RT + v2);
	M16[13] =  ((degrees + 4.0)*RT + u2 + v2)*M16[4];
	M16[14] =  degPlus2*rho*RT*RT + p*((degrees + 5.0)*u2 + v2) + rho*u2*(u2 + v2);
	M16[15] =  degPlus2*rho*RT*RT + p*((degrees + 5.0)*v2 + u2) + rho*v2*(u2 + v2);
}
double vectorDot(const double *const first,const double *const second)
{
	double sum = 0;
	for(int j = 0;j < DV_Qv;++j)
		sum += first[j]*second[j];
	return sum;
}
double vectorDot(const double *const first,const double *const second,const double *const third)
{
	double sum = 0;
	for(int j = 0;j < DV_Qv;++j)
		sum += first[j]*second[j]*third[j];
	return sum;
}
double vectorDot(const double *const first)
{
	double sum = 0;
	for(int j = 0;j < DV_Qv;++j)
		sum += first[j];
	return sum;
}
void Update_phi_Eq(Cell_2D &cell)
{
	double M16[DV_Qv];
	Macro16(M16,cell.Rho,cell.U,cell.V,cell.T,cell.p);
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		cell.fEq[i][j] = vectorDot(iArkC[j],M16);
	}
}
void Update_phi_Eqh(Face_2D &face)
{
	double M16[DV_Qv];
	Macro16(M16,face.Rho_h,face.U_h,face.V_h,face.T_h,face.p_h);
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		face.fEqh[i][j] = vectorDot(iArkC[j],M16);
	}
}
void Update_MacroVar(Cell_2D& cell)
{
	cell.Rho  = vectorDot(cell.fT[0]);
	cell.U    = vectorDot(cell.fT[0],xi_u)/cell.Rho;
	cell.V    = vectorDot(cell.fT[0],xi_v)/cell.Rho;
	cell.T    = vectorDot(cell.fT[0],xi2Eta2)/cell.Rho;
	cell.T = (cell.T - cell.U*cell.U - cell.V*cell.V)/degrees;
	cell.p = cell.Rho*R0*cell.T;
	cell.Mu = Mu0*pow(cell.T/T0,Omega0);
	cell.Factor();
}

void Update_MacroVar_h(Face_2D& face)
{
	face.Rho_h  = vectorDot(face.fBh[0]);
	face.U_h    = vectorDot(face.fBh[0],xi_u)/face.Rho_h;
	face.V_h    = vectorDot(face.fBh[0],xi_v)/face.Rho_h;
	face.T_h    = vectorDot(face.fBh[0],xi2Eta2)/face.Rho_h;
	face.T_h = (face.T_h - face.U_h*face.U_h - face.V_h*face.V_h)/degrees;
	face.p_h = face.Rho_h*R0*face.T_h;
	face.Mu_h = Mu0*pow(face.T_h/T0,Omega0);
	face.Factor();
}
void IntegralShearStress()
{
	for(int n = 0;n < Cells;++n)
	{
		double fNeq[DV_Qv];
		Cell_2D &cell = CellArray[n];
		double ashearTau = cell.Tau/(cell.Tau + h);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			fNeq[j] = cell.fT[i][j] - cell.fEq[i][j];
		}
		cell.shearTau[0][0] = ashearTau*vectorDot(fNeq,xi_u,xi_u);
		cell.shearTau[0][1] = ashearTau*vectorDot(fNeq,xi_u,xi_v);
		cell.shearTau[1][0] = cell.shearTau[0][1];
		cell.shearTau[1][1] = ashearTau*vectorDot(fNeq,xi_v,xi_v);
	}
}
void UW_Interior_phi_Bh_Limiter(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
{
	double dx = face.xf - h*xi_u[j] - ptr_C->xc;
	double dy = face.yf - h*xi_v[j] - ptr_C->yc;
	VenkatakrishnanFluxLimiter(*ptr_C,i,j);
	face.fBh[i][j] = ptr_C->fBP[i][j] + ptr_C->fBPLimiter*(dx*ptr_C->fBP_x[i][j] + dy*ptr_C->fBP_y[i][j]);
}
void UW_Interior_phi_Bh(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
{
	double dx = face.xf - h*xi_u[j] - ptr_C->xc;
	double dy = face.yf - h*xi_v[j] - ptr_C->yc;
	face.fBh[i][j] = ptr_C->fBP[i][j] + (dx*ptr_C->fBP_x[i][j] + dy*ptr_C->fBP_y[i][j]);
}
void CD_Interior_phi_Bh(Face_2D &face,int i,int j)
{
	double _dx = face.lhsCell->xc - face.rhsCell->xc;
	double _dy = face.lhsCell->yc - face.rhsCell->yc;
	SetZero(_dx);
	SetZero(_dy);
	double _fBP_xF,_fBP_yF;
	{
		if(0.0 == _dx)
		{
			 _fBP_xF =  0.5*(face.lhsCell->fBP_x[i][j] + face.rhsCell->fBP_x[i][j]);
			 _fBP_yF = (face.lhsCell->fBP[i][j] - face.rhsCell->fBP[i][j])/_dy;
		}
		else if(0.0 == _dy)
		{
			_fBP_yF = 0.5*(face.lhsCell->fBP_y[i][j] + face.rhsCell->fBP_y[i][j]);
			_fBP_xF = (face.lhsCell->fBP[i][j] - face.rhsCell->fBP[i][j])/_dx;
		}
		else
		{
			_PRINT_ERROR_MSG_FLIP
			getchar();
		}
		face.fBh[i][j] = 0.5*(face.lhsCell->fBP[i][j] + face.rhsCell->fBP[i][j])
				 - h*(_fBP_xF*xi_u[j] + _fBP_yF*xi_v[j]);
	}	
}