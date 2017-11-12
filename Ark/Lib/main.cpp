#include <iostream>
#include <cmath>
#include <iomanip>
using std::cout;
using std::setprecision;
#include "NewtonCotes.h"

int const 

DV_Qu = 101, 

DV_Qv = 101;

double const

Rho0 = 1,

U0 = 0,

V0 = 0,

T0 = 0.5,

R0 = 1,

Lambda = 1.0/(2*R0*T0),

nK = 2,

agEq = (nK + 3.0 - 2.0)/2.0,

Gamma = (nK + 5.0)/(nK + 3.0);

double const

PI = 3.14159265358,

MaSpan = 10.0,

MaxU = 2*MaSpan*sqrt(R0*T0),

UpperLimitQu = MaSpan*sqrt(2.0*R0*T0),

LowerLimitQu = -UpperLimitQu,//0,

UpperLimitQv = MaSpan*sqrt(2.0*R0*T0),

LowerLimitQv = -UpperLimitQv,

DeltaQu = (UpperLimitQu - LowerLimitQu)/(DV_Qu - 1.0),

DeltaQv = (UpperLimitQv - LowerLimitQv)/(DV_Qv - 1.0);

double **fEq,**gEq;

double Rho,U,V,p,E,T,xi_u[DV_Qu],xi_v[DV_Qv];

void DiscreteVelocityAssign();

void AllocateARK(double** &f,int const Qu,int const Qv);

void DeallocateARK(double** &f,int const Qu,int const Qv);

void Update_phi_Eq();

void Update_MacroVar();

int main()
{
	DiscreteVelocityAssign();
//
	AllocateARK(fEq,DV_Qu,DV_Qv);
	AllocateARK(gEq,DV_Qu,DV_Qv);

	 Update_phi_Eq();
	 Update_MacroVar();
	// double f[DV_Qu];
	// cout <<DeltaQu<<'\n';
	// for(int i = 0;i < DV_Qu;++i)
	// 	f[i] = sqrt(xi_u[i]);
	// Rho = CompositeCotes(f,DV_Qu,DeltaQu);
	// cout<<setprecision(16)<<Rho<<'\n';

	DeallocateARK(fEq,DV_Qu,DV_Qv);
	DeallocateARK(gEq,DV_Qu,DV_Qv);
}

//-------------------------------------Preprocess-----------------------------
void AllocateARK(double** &f,int const Qu,int const Qv)
{
	f = new double* [Qu];
	for(int i = 0;i < Qu;++i)
		f[i] = new double[Qv];
}
void DeallocateARK(double** &f,int const Qu,int const Qv)
{
	for(int i = 0;i < Qu;++i)
		delete[] f[i];
	delete[] f;
}
void DiscreteVelocityAssign()
{
	for(int i = 0;i < DV_Qu;++i)
	{
		xi_u[i] = LowerLimitQu + i*DeltaQu;
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		xi_v[j] = LowerLimitQv + j*DeltaQv;
	}
}
void Update_phi_Eq()
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - U0, cy = xi_v[j] - V0;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-Lambda*c2)/PI;
		fEq[i][j] = Lambda*Rho0*aEq;
		gEq[i][j] = agEq*Rho0*aEq;
		cout <<fEq[i][j]<<'\n';
	}
}
void Update_MacroVar()
{
	double Int_fEq_du[DV_Qu],Int_fEq_dv[DV_Qv],Int_gEq_du[DV_Qu];
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_fEq_du[i] = CompositeCotes(fEq[i],DV_Qv,DeltaQv);
		Int_gEq_du[i] = CompositeCotes(gEq[i],DV_Qv,DeltaQv);

		cout <<i<<"  "<<Int_fEq_du[i]<<'\n';
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		Int_fEq_dv[j] = CompositeCotes(fEq,DV_Qu,j,DeltaQu);
	}
	Rho = CompositeCotes(Int_fEq_du,DV_Qu,DeltaQu);
	U   = CompositeCotes(Int_fEq_du,DV_Qu,xi_u,DeltaQu)/Rho;
	V   = CompositeCotes(Int_fEq_dv,DV_Qv,xi_v,DeltaQv)/Rho;
	E   = (CompositeCotes(Int_gEq_du,DV_Qu,DeltaQu) + CompositeCotes(Int_fEq_du,DV_Qu,xi_u,xi_u,DeltaQu)
				+ CompositeCotes(Int_fEq_dv,DV_Qv,xi_v,xi_v,DeltaQv))*0.5;
	E   -= 0.5*Rho*(U*U + V*V);
	p = E*(Gamma - 1.0);
	T = p/(Rho*R0);
	cout <<"Rho : "<<Rho<<'\n'
		 <<U<<'\n'
		 <<V<<'\n'
		 <<T<<'\n';
}
