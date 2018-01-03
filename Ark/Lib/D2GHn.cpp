#include <fstream>
#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

using std::ifstream;

extern double * const xi_u = new double[DV_Qu];

extern double * const xi_v = new double[DV_Qv];

double wGH[DV_Qu];

double const 

UpperLimitQu = MaSpan*sqrt(2.0*R0*T0),

LowerLimitQu = -UpperLimitQu,

UpperLimitQv = MaSpan*sqrt(2.0*R0*T0),

LowerLimitQv = -UpperLimitQv;

extern void AllocateARK(double** &f,int const Qu,int const Qv);

extern void DeallocateARK(double** &f,int const Qu,int const Qv);

void DiscreteVelocityAssign()
{
	int const Lines = 128;
	char cc[Lines] = {'0'};
	ifstream InFile_Xis("./constant/Xis");
	ifstream InFile_weights("./constant/weights");
	while('(' != cc[0])
	{
		InFile_Xis.getline(cc,Lines);
	}
	cc[0] = 'o';
	while('(' != cc[0])
	{
		InFile_weights.getline(cc,Lines);
	}
	for(int i = 0;i < DV_Qu;++i)
	{
		InFile_Xis >> xi_u[i];
		InFile_weights >> wGH[i];
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		xi_v[j] = xi_u[j];
	}
}

void setXiDotdS()
{
	for(int n = 0;n < Faces;++n)
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		FaceArray[n].xi_n_dS[i][j] = FaceArray[n].Area*(xi_u[i]*FaceArray[n].Vx + xi_v[j]*FaceArray[n].Vy);
	}
}
void Update_phi_Eq(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - cell.U, cy = xi_v[j] - cell.V;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-cell.Lambda*c2)/PI;
#ifndef _ARK_ISOTHERMAL_FLIP
		double CQ = (1.0-Pr)*cell.Lambda*cell.Lambda*(cx*cell.qx+cy*cell.qy);
		cell.fEq[i][j] = (1.6*(cell.Lambda*c2-2.0)*CQ + cell.Rho)*cell.Lambda*aEq;
		cell.gEq[i][j] = (0.8*((cell.Lambda*c2-1.0)*(nK+1.0)-nK)*CQ + agEq*cell.Rho)*aEq;
#else
		cell.fEq[i][j] = cell.Lambda*cell.Rho*aEq;
#endif
	}
}

void Update_phi_Eqh(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - face.U_h, cy = xi_v[j] - face.V_h;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-face.Lambda_h*c2)/PI;
#ifndef _ARK_ISOTHERMAL_FLIP
		double CQ = (1.0-Pr)*face.Lambda_h*face.Lambda_h*(cx*face.qx_h + cy*face.qy_h);
		face.fEqh[i][j] = (1.6*(face.Lambda_h*c2-2.0)*CQ + face.Rho_h)*face.Lambda_h*aEq;
		face.gEqh[i][j] = (0.8*((face.Lambda_h*c2-1.0)*(nK+1.0)-nK)*CQ + agEq*face.Rho_h)*aEq;
#else
		face.fEqh[i][j] = face.Rho_h*face.Lambda_h*aEq;
#endif
	}
}

void Update_MacroVar_h(Face_2D& face)
{
	double Int_fBh_du[DV_Qu],Int_fBh_dv[DV_Qv],Int_gBh_du[DV_Qu];
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_fBh_du[i] = IntegralGH(face.fBh[i],DV_Qv,wGH);
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		Int_fBh_dv[j] = IntegralGH(face.fBh,DV_Qu,j,wGH);
	}
	face.Rho_h = IntegralGH(Int_fBh_du,DV_Qu,wGH);
	face.U_h   = IntegralGH(Int_fBh_du,DV_Qu,xi_u,wGH)/face.Rho_h;
	face.V_h   = IntegralGH(Int_fBh_dv,DV_Qv,xi_v,wGH)/face.Rho_h;
#ifndef _ARK_ISOTHERMAL_FLIP
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_gBh_du[i] = IntegralGH(face.gBh[i],DV_Qv,wGH);
	}
	face.E_h   = (IntegralGH(Int_gBh_du,DV_Qu,wGH) + IntegralGH(Int_fBh_du,DV_Qu,xi_u,xi_u,wGH)
				+ IntegralGH(Int_fBh_dv,DV_Qv,xi_v,xi_v,wGH))*0.5;
	face.E_h  -= 0.5*face.Rho_h*(face.U_h*face.U_h + face.V_h*face.V_h);
	face.p_h = face.E_h*(Gamma - 1.0);
	face.T_h = face.p_h/(face.Rho_h*R0);
	face.Lambda_h = 0.5/(R0*face.T_h);
	face.Mu_h = Mu0*pow(face.T_h/T0,Omega0);
	face.Factor();
//
	double Int_qxh_du[DV_Qu],Int_qyh_du[DV_Qu],Cx[DV_Qu],Cy[DV_Qv];
	double **Cx_fBhC2_gBh,**Cy_fBhC2_gBh,fBhC2_gBh;
	AllocateARK(Cx_fBhC2_gBh,DV_Qu,DV_Qv);
	AllocateARK(Cy_fBhC2_gBh,DV_Qu,DV_Qv);
	for(int i = 0;i < DV_Qu;++i)
		Cx[i] = xi_u[i] - face.U_h;
	for(int j = 0;j < DV_Qv;++j)
		Cy[j] = xi_v[j] - face.V_h;
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		fBhC2_gBh = (Cx[i]*Cx[i] + Cy[j]*Cy[j])*face.fBh[i][j] + face.gBh[i][j];
		Cx_fBhC2_gBh[i][j] = fBhC2_gBh*Cx[i];
		Cy_fBhC2_gBh[i][j] = fBhC2_gBh*Cy[j];
	}
	for(int i = 0;i < DV_Qu;++i)
	{
		Int_qxh_du[i] = IntegralGH(Cx_fBhC2_gBh[i],DV_Qv,wGH);
		Int_qyh_du[i] = IntegralGH(Cy_fBhC2_gBh[i],DV_Qv,wGH);
	}
	face.qx_h = face.aQh*IntegralGH(Int_qxh_du,DV_Qu,wGH);
	face.qy_h = face.aQh*IntegralGH(Int_qyh_du,DV_Qu,wGH);
	DeallocateARK(Cx_fBhC2_gBh,DV_Qu,DV_Qv);
	DeallocateARK(Cy_fBhC2_gBh,DV_Qu,DV_Qv);
#endif
}

void Update_MacroVar(Cell_2D& cell)
{
	double Int_fT_du[DV_Qu],Int_fT_dv[DV_Qv],Int_gT_du[DV_Qu];
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_fT_du[i] = IntegralGH(cell.fT[i],DV_Qv,wGH);
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		Int_fT_dv[j] = IntegralGH(cell.fT,DV_Qu,j,wGH);
	}
	cell.Rho = IntegralGH(Int_fT_du,DV_Qu,wGH);
	cell.U   = IntegralGH(Int_fT_du,DV_Qu,xi_u,wGH)/cell.Rho;
	cell.V   = IntegralGH(Int_fT_dv,DV_Qv,xi_v,wGH)/cell.Rho;
#ifndef _ARK_ISOTHERMAL_FLIP
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_gT_du[i] = IntegralGH(cell.gT[i],DV_Qv,wGH);
	}
	cell.E   = (IntegralGH(Int_gT_du,DV_Qu,wGH) + IntegralGH(Int_fT_du,DV_Qu,xi_u,xi_u,wGH)
				+ IntegralGH(Int_fT_dv,DV_Qv,xi_v,xi_v,wGH))*0.5;
	cell.E  -= 0.5*cell.Rho*(cell.U*cell.U + cell.V*cell.V);
	cell.p = cell.E*(Gamma - 1.0);
	cell.T = cell.p/(cell.Rho*R0);
	cell.Lambda = 0.5/(R0*cell.T);
	cell.Mu = Mu0*pow(cell.T/T0,Omega0);
	cell.Factor();
//
	double Int_qx_du[DV_Qu],Int_qy_du[DV_Qu],Cx[DV_Qu],Cy[DV_Qv];
	double **Cx_fTC2_gT,**Cy_fTC2_gT,fTC2_gT;
	AllocateARK(Cx_fTC2_gT,DV_Qu,DV_Qv);
	AllocateARK(Cy_fTC2_gT,DV_Qu,DV_Qv);
	for(int i = 0;i < DV_Qu;++i)
		Cx[i] = xi_u[i] - cell.U;
	for(int j = 0;j < DV_Qv;++j)
		Cy[j] = xi_v[j] - cell.V;
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		fTC2_gT = (Cx[i]*Cx[i] + Cy[j]*Cy[j])*cell.fT[i][j] + cell.gT[i][j];
		Cx_fTC2_gT[i][j] = fTC2_gT*Cx[i];
		Cy_fTC2_gT[i][j] = fTC2_gT*Cy[j];
	}
	for(int i = 0;i < DV_Qu;++i)
	{
		Int_qx_du[i] = IntegralGH(Cx_fTC2_gT[i],DV_Qv,wGH);
		Int_qy_du[i] = IntegralGH(Cy_fTC2_gT[i],DV_Qv,wGH);
	}
	cell.qx = cell.aQ*IntegralGH(Int_qx_du,DV_Qu,wGH);
	cell.qy = cell.aQ*IntegralGH(Int_qy_du,DV_Qu,wGH);
	DeallocateARK(Cx_fTC2_gT,DV_Qu,DV_Qv);
	DeallocateARK(Cy_fTC2_gT,DV_Qu,DV_Qv);
#endif
}
void UW_Interior_phi_Bh_Limiter(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
{
	double dx = face.xf - h*xi_u[i] - ptr_C->xc;
	double dy = face.yf - h*xi_v[j] - ptr_C->yc;
	VenkatakrishnanFluxLimiter(*ptr_C,i,j);
	face.fBh[i][j] = ptr_C->fBP[i][j] + ptr_C->fBPLimiter*(dx*ptr_C->fBP_x[i][j] + dy*ptr_C->fBP_y[i][j]);
//isothermal flip
	#ifndef _ARK_ISOTHERMAL_FLIP
	face.gBh[i][j] = ptr_C->gBP[i][j] + ptr_C->gBPLimiter*(dx*ptr_C->gBP_x[i][j] + dy*ptr_C->gBP_y[i][j]);
	#endif
}
void UW_Interior_phi_Bh(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
{
	double dx = face.xf - h*xi_u[i] - ptr_C->xc;
	double dy = face.yf - h*xi_v[j] - ptr_C->yc;
	face.fBh[i][j] = ptr_C->fBP[i][j] + (dx*ptr_C->fBP_x[i][j] + dy*ptr_C->fBP_y[i][j]);
//isothermal flip
	#ifndef _ARK_ISOTHERMAL_FLIP
	face.gBh[i][j] = ptr_C->gBP[i][j] + (dx*ptr_C->gBP_x[i][j] + dy*ptr_C->gBP_y[i][j]);
	#endif
}
void CD_Interior_phi_Bh(Face_2D &face,int i,int j)
{
	double _dx = face.lhsCell->xc - face.rhsCell->xc;
	double _dy = face.lhsCell->yc - face.rhsCell->yc;
	SetZero(_dx);
	SetZero(_dy);
	double _fBP_xF,_fBP_yF,_gBP_xF,_gBP_yF;
	{
		if(0.0 == _dx)
		{
			 _fBP_xF =  0.5*(face.lhsCell->fBP_x[i][j] + face.rhsCell->fBP_x[i][j]);
			 _fBP_yF = (face.lhsCell->fBP[i][j] - face.rhsCell->fBP[i][j])/_dy;
			 #ifndef _ARK_ISOTHERMAL_FLIP
			 _gBP_xF =  0.5*(face.lhsCell->gBP_x[i][j] + face.rhsCell->gBP_x[i][j]);
			 _gBP_yF = (face.lhsCell->gBP[i][j] - face.rhsCell->gBP[i][j])/_dy;
			 #endif
		}
		else if(0.0 == _dy)
		{
			_fBP_yF = 0.5*(face.lhsCell->fBP_y[i][j] + face.rhsCell->fBP_y[i][j]);
			_fBP_xF = (face.lhsCell->fBP[i][j] - face.rhsCell->fBP[i][j])/_dx;
			#ifndef _ARK_ISOTHERMAL_FLIP
			_gBP_yF = 0.5*(face.lhsCell->gBP_y[i][j] + face.rhsCell->gBP_y[i][j]);
			_gBP_xF = (face.lhsCell->gBP[i][j] - face.rhsCell->gBP[i][j])/_dx;
			#endif
		}
		else
		{
		}
		face.fBh[i][j] = 0.5*(face.lhsCell->fBP[i][j] + face.rhsCell->fBP[i][j])
				 - h*(_fBP_xF*xi_u[i] + _fBP_yF*xi_v[j]);
		#ifndef _ARK_ISOTHERMAL_FLIP
		face.gBh[i][j] = 0.5*(face.lhsCell->gBP[i][j] + face.rhsCell->gBP[i][j])
				 - h*(_gBP_xF*xi_u[i] + _gBP_yF*xi_v[j]);
		#endif
	}	
}
void IntegralShearStress()
{
	for(int n = 0;n < Cells;++n)
	{
		Cell_2D& cell = CellArray[n];
		double fTMinusfEq,Int_sTau_xx_du[DV_Qu],Int_sTau_xy_du[DV_Qu],Int_sTau_yy_du[DV_Qu];
		double **sTau_xx,**sTau_xy,**sTau_yy;//xx,xy,yx,yy
		double ashearTau = cell.Tau/(cell.Tau + h);
		AllocateARK(sTau_xx,DV_Qu,DV_Qv);
		AllocateARK(sTau_xy,DV_Qu,DV_Qv);
		AllocateARK(sTau_yy,DV_Qu,DV_Qv);
		Update_phi_Eq(cell);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			fTMinusfEq = cell.fT[i][j] - cell.fEq[i][j];
			sTau_xx[i][j] = xi_u[i]*xi_u[i]*fTMinusfEq;
			sTau_xy[i][j] = xi_u[i]*xi_v[j]*fTMinusfEq;
			sTau_yy[i][j] = xi_v[j]*xi_v[j]*fTMinusfEq;
		}
		for(int i = 0;i < DV_Qu;++i)
		{
			Int_sTau_xx_du[i] = IntegralGH(sTau_xx[i],DV_Qv,wGH);
			Int_sTau_xy_du[i] = IntegralGH(sTau_xy[i],DV_Qv,wGH);
			Int_sTau_yy_du[i] = IntegralGH(sTau_yy[i],DV_Qv,wGH);
		}
		cell.shearTau[0][0] = ashearTau*IntegralGH(Int_sTau_xx_du,DV_Qu,wGH);
		cell.shearTau[0][1] = ashearTau*IntegralGH(Int_sTau_xy_du,DV_Qu,wGH);
		cell.shearTau[1][0] = cell.shearTau[0][1];
		cell.shearTau[1][1] = ashearTau*IntegralGH(Int_sTau_yy_du,DV_Qu,wGH);
		DeallocateARK(sTau_xx,DV_Qu,DV_Qv);
		DeallocateARK(sTau_xy,DV_Qu,DV_Qv);
		DeallocateARK(sTau_yy,DV_Qu,DV_Qv);
	}
}
//---------------------------------Wall_3_Boundary------------------------
// void Wall_3_Boundary(Face_2D &face)
// {
// 	double sumOut = 0,sumIn = 0;
// //-----------------------Real Flux Outward---------------------
// 	face.Rho_h = face.lhsCell->Rho;
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		if(face.xi_n_dS[i][j] > 0)
// 		{
// 			UW_Interior_phi_Bh(face,face.lhsCell,i,j);
// 			Update_phi_Eqh(face,i,j);
// 			Update_phi_h(face,i,j);
// 			Update_phiFlux_h(face,i,j);
// 			sumOut += wGH[i]*wGH[j]*face.fh[i][j];
// 		}
// 	}
// //---------------------Spurious Flux Inward-----------------------
// 	face.Rho_h = 1.0;
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		if(face.xi_n_dS[i][j] < 0)
// 		{
// 			Update_phi_Eqh(face,i,j);
// 			face.fh[i][j] = face.fEqh[i][j];
// 			//isothermal flip
// 			#ifndef _ARK_ISOTHERMAL_FLIP
// 			face.gh[i][j] = face.gEqh[i][j];
// 			#endif
// 			Update_phiFlux_h(face,i,j);
// 			sumIn += wGH[i]*wGH[j]*face.fh[i][j];
// 		}
// 	}
// 	face.Rho_h = -sumOut/sumIn;
// //------------------------------------------------
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		if(face.xi_n_dS[i][j] < 0)
// 		{
// 			Update_phi_Eqh(face,i,j);
// 			face.fh[i][j] = face.fEqh[i][j];
// 			//isothermal flip
// 			#ifndef _ARK_ISOTHERMAL_FLIP
// 			face.gh[i][j] = face.gEqh[i][j];
// 			#endif
// 			Update_phiFlux_h(face,i,j);
// 		}
// 	}
// }