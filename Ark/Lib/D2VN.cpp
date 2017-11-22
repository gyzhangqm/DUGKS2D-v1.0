#include "DUGKSDeclaration.h"
#include "NewtonCotes.h"

extern double * const xi_u = new double[DV_Qu];

extern double * const xi_v = new double[DV_Qv];

double const 

agEq = (nK + 3.0 - 2.0)/2.0,

UpperLimitQu = MaSpan*sqrt(2.0*R0*T0),

LowerLimitQu = -UpperLimitQu,

UpperLimitQv = MaSpan*sqrt(2.0*R0*T0),

LowerLimitQv = -UpperLimitQv,

DeltaQu = (UpperLimitQu - LowerLimitQu)/(DV_Qu - 1.0),

DeltaQv = (UpperLimitQv - LowerLimitQv)/(DV_Qv - 1.0);

extern void AllocateARK(double** &f,int const Qu,int const Qv);

extern void DeallocateARK(double** &f,int const Qu,int const Qv);

//-------------------------------------Preprocess-----------------------------
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
//
void setXiDotdS()
{
	for(int n = 0;n < Faces;++n)
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		FaceArray[n].xi_n_dS[i][j] = FaceArray[n].Area*(xi_u[i]*FaceArray[n].Vx + xi_v[j]*FaceArray[n].Vy);
	}
}
//-------------------------------------------------------------------------

//-------------------------------------Shakhov Equilibrium-----------------
/*void Update_phi_Eq(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - cell.U, cy = xi_v[j] - cell.V;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-cell.Lambda*c2)/PI;
		cell.fEq[i][j] = cell.Lambda*cell.Rho*aEq;
		cell.gEq[i][j] = agEq*cell.Rho*aEq;
	}
}*/
void Update_phi_Eq(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - cell.U, cy = xi_v[j] - cell.V;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-cell.Lambda*c2)/PI;
		double CQ = (1.0-Pr)*cell.Lambda*cell.Lambda*(cx*cell.qx+cy*cell.qy);
		cell.fEq[i][j] = (1.6*(cell.Lambda*c2-2.0)*CQ + cell.Rho)*cell.Lambda*aEq;
		cell.gEq[i][j] = (0.8*((cell.Lambda*c2-1.0)*(nK+1.0)-nK)*CQ + agEq*cell.Rho)*aEq;
	}
}
/*void Update_phi_Eqh(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - face.U_h, cy = xi_v[j] - face.V_h;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-face.Lambda_h*c2)/PI;
		face.fEqh[i][j] = face.Rho_h*face.Lambda_h*aEq;
		face.gEqh[i][j] = agEq*face.Rho_h*aEq;
	}
}*/
void Update_phi_Eqh(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		double cx = xi_u[i] - face.U_h, cy = xi_v[j] - face.V_h;
		double c2 = cx*cx + cy*cy;
		double aEq = exp(-face.Lambda_h*c2)/PI;
		double CQ = (1.0-Pr)*face.Lambda_h*face.Lambda_h*(cx*face.qx_h + cy*face.qy_h);
		face.fEqh[i][j] = (1.6*(face.Lambda_h*c2-2.0)*CQ + face.Rho_h)*face.Lambda_h*aEq;
		face.gEqh[i][j] = (0.8*((face.Lambda_h*c2-1.0)*(nK+1.0)-nK)*CQ + agEq*face.Rho_h)*aEq;
	}
}
//-----------------------------------Newton-Cotes---------------------------------
void Update_MacroVar(Cell_2D& cell)
{
	double Int_fT_du[DV_Qu],Int_fT_dv[DV_Qv],Int_gT_du[DV_Qu];
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_fT_du[i] = CompositeCotes(cell.fT[i],DV_Qv,DeltaQv);
		Int_gT_du[i] = CompositeCotes(cell.gT[i],DV_Qv,DeltaQv);
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		Int_fT_dv[j] = CompositeCotes(cell.fT,DV_Qu,j,DeltaQu);
	}
	cell.Rho = CompositeCotes(Int_fT_du,DV_Qu,DeltaQu);
	cell.U   = CompositeCotes(Int_fT_du,DV_Qu,xi_u,DeltaQu)/cell.Rho;
	cell.V   = CompositeCotes(Int_fT_dv,DV_Qv,xi_v,DeltaQv)/cell.Rho;
	cell.E   = (CompositeCotes(Int_gT_du,DV_Qu,DeltaQu) + CompositeCotes(Int_fT_du,DV_Qu,xi_u,xi_u,DeltaQu)
				+ CompositeCotes(Int_fT_dv,DV_Qv,xi_v,xi_v,DeltaQv))*0.5;
	cell.E   -= 0.5*cell.Rho*(cell.U*cell.U + cell.V*cell.V);
	cell.p = cell.E*(Gamma - 1.0);
	cell.T = cell.p/(cell.Rho*R0);
	cell.Lambda = 0.5/(R0*cell.T);
	cell.Mu = Mu0*pow(cell.T/T0,Omega0);
	cell.Factor();
//heat flux
	double Int_qx_du[DV_Qu],Int_qy_du[DV_Qu],Cx[DV_Qu],Cy[DV_Qv];
	double **Cx_fTC2_gT,**Cy_fTC2_gT,fTC2_gT;
//
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
		Int_qx_du[i] = CompositeCotes(Cx_fTC2_gT[i],DV_Qv,DeltaQv);
		Int_qy_du[i] = CompositeCotes(Cy_fTC2_gT[i],DV_Qv,DeltaQv);
	}
	cell.qx = cell.aQ*CompositeCotes(Int_qx_du,DV_Qu,DeltaQu);
	cell.qy = cell.aQ*CompositeCotes(Int_qy_du,DV_Qu,DeltaQu);
	DeallocateARK(Cx_fTC2_gT,DV_Qu,DV_Qv);
	DeallocateARK(Cy_fTC2_gT,DV_Qu,DV_Qv);
}
void Update_MacroVar_h(Face_2D& face)
{
	double Int_fBh_du[DV_Qu],Int_fBh_dv[DV_Qv],Int_gBh_du[DV_Qu];
	for(int i = 0;i < DV_Qu;++i)		
	{
		Int_fBh_du[i] = CompositeCotes(face.fBh[i],DV_Qv,DeltaQv);
		Int_gBh_du[i] = CompositeCotes(face.gBh[i],DV_Qv,DeltaQv);
	}
	for(int j = 0;j < DV_Qv;++j)
	{
		Int_fBh_dv[j] = CompositeCotes(face.fBh,DV_Qu,j,DeltaQu);
	}
	face.Rho_h = CompositeCotes(Int_fBh_du,DV_Qu,DeltaQu);
	face.U_h   = CompositeCotes(Int_fBh_du,DV_Qu,xi_u,DeltaQu)/face.Rho_h;
	face.V_h   = CompositeCotes(Int_fBh_dv,DV_Qv,xi_v,DeltaQv)/face.Rho_h;
	face.E_h   = 0.5*(CompositeCotes(Int_fBh_du,DV_Qu,xi_u,xi_u,DeltaQu)
				 	+ CompositeCotes(Int_fBh_dv,DV_Qv,xi_v,xi_v,DeltaQv)
				 	+ CompositeCotes(Int_gBh_du,DV_Qu,DeltaQu));
	face.E_h  -= 0.5*face.Rho_h*(face.U_h*face.U_h + face.V_h*face.V_h);
	face.p_h = face.E_h*(Gamma - 1.0);
	face.T_h = face.p_h/(face.Rho_h*R0);
	face.Lambda_h = 0.5/(R0*face.T_h);
	face.Mu_h = Mu0*pow(face.T_h/T0,Omega0);
	face.Factor();
//heat flux
	double Int_qxh_du[DV_Qu],Int_qyh_du[DV_Qu],Cx[DV_Qu],Cy[DV_Qv];
	double **Cx_fBhC2_gBh,**Cy_fBhC2_gBh,fBhC2_gBh;
//
	AllocateARK(Cx_fBhC2_gBh,DV_Qu,DV_Qv);
	AllocateARK(Cy_fBhC2_gBh,DV_Qu,DV_Qv);
//
	for(int i = 0;i < DV_Qu;++i)
		Cx[i] = xi_u[i] - face.U_h;
	for(int j = 0;j < DV_Qv;++j)
		Cy[j] = xi_v[j] - face.V_h;
//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		fBhC2_gBh = (Cx[i]*Cx[i] + Cy[j]*Cy[j])*face.fBh[i][j] + face.gBh[i][j];
		Cx_fBhC2_gBh[i][j] = Cx[i]*fBhC2_gBh;
		Cy_fBhC2_gBh[i][j] = Cy[j]*fBhC2_gBh;
	}
	for(int i = 0;i < DV_Qu;++i)
	{
		Int_qxh_du[i] = CompositeCotes(Cx_fBhC2_gBh[i],DV_Qv,DeltaQv);
		Int_qyh_du[i] = CompositeCotes(Cy_fBhC2_gBh[i],DV_Qv,DeltaQv);
	}
	face.qx_h = face.aQh*CompositeCotes(Int_qxh_du,DV_Qu,DeltaQu);
	face.qy_h = face.aQh*CompositeCotes(Int_qyh_du,DV_Qu,DeltaQu);
//
	DeallocateARK(Cx_fBhC2_gBh,DV_Qu,DV_Qv);
	DeallocateARK(Cy_fBhC2_gBh,DV_Qu,DV_Qv);
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
			Int_sTau_xx_du[i] = CompositeCotes(sTau_xx[i],DV_Qv,DeltaQv);
			Int_sTau_xy_du[i] = CompositeCotes(sTau_xy[i],DV_Qv,DeltaQv);
			Int_sTau_yy_du[i] = CompositeCotes(sTau_yy[i],DV_Qv,DeltaQv);
		}
		cell.shearTau[0][0] = ashearTau*CompositeCotes(Int_sTau_xx_du,DV_Qu,DeltaQu);
		cell.shearTau[0][1] = ashearTau*CompositeCotes(Int_sTau_xy_du,DV_Qu,DeltaQu);
		cell.shearTau[1][0] = cell.shearTau[0][1];
		cell.shearTau[1][1] = ashearTau*CompositeCotes(Int_sTau_yy_du,DV_Qu,DeltaQu);
		DeallocateARK(sTau_xx,DV_Qu,DV_Qv);
		DeallocateARK(sTau_xy,DV_Qu,DV_Qv);
		DeallocateARK(sTau_yy,DV_Qu,DV_Qv);
	}
}