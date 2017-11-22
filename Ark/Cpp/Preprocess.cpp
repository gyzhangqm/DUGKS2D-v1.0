#include <iostream>
#include "DUGKSDeclaration.h"
//
using std::cout;
using std::endl;

//--------------------------------------------------------------------------------------
//------------------------------------------Initialization------------------------------
void UniformFlow()
{
		for(int n = 0;n < Cells;++n)
		{
			CellArray[n].Rho = Rho_Outlet;
			CellArray[n].U = U_Outlet;
			CellArray[n].V = V_Outlet;
			CellArray[n].T = T_Outlet;
			CellArray[n].p = Rho_Outlet*R0*T_Outlet;
			CellArray[n].Lambda = 0.5/(T_Outlet*R0);
			CellArray[n].qx = 0;
			CellArray[n].qy = 0;
			CellArray[n].Mu = Mu0;
			CellArray[n].Factor();
			Update_phi_Eq(CellArray[n]);
			for(int i = 0;i < DV_Qu;++i)
			for(int j = 0;j < DV_Qv;++j)
			{
				CellArray[n].fT[i][j]  = CellArray[n].fEq[i][j];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
				CellArray[n].gT[i][j]  = CellArray[n].gEq[i][j];
#endif
			}
		}
		for(int n = 0;n < Faces;++n)
		{
			FaceArray[n].Rho_h = Rho_Outlet;
			FaceArray[n].U_h = U_Outlet;
			FaceArray[n].V_h = V_Outlet;
			FaceArray[n].T_h = T_Outlet;
			FaceArray[n].p_h = Rho_Outlet*R0*T_Outlet;
			FaceArray[n].Lambda_h = 0.5/(T_Outlet*R0);
			FaceArray[n].Mu_h = Mu0;
			FaceArray[n].Factor();
			for(int i = 0;i < DV_Qu;++i)
			for(int j = 0;j < DV_Qv;++j)
			{
				FaceArray[n].fEqh[i][j] = FaceArray[n].lhsCell->fEq[i][j];
				FaceArray[n].fBh[i][j] = FaceArray[n].fEqh[i][j];
				FaceArray[n].fh[i][j] = FaceArray[n].fEqh[i][j];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
				FaceArray[n].gEqh[i][j] = FaceArray[n].lhsCell->gEq[i][j];
				FaceArray[n].gBh[i][j] = FaceArray[n].gEqh[i][j];
				FaceArray[n].gh[i][j] = FaceArray[n].gEqh[i][j];
#endif
			}
		}
		#ifdef _P_INLET_4_BCS_FLIP
		for(int k = 0;k < P_InletFaceNum;++k)
		{
			for(int i = 0;i < DV_Qu;++i)
			for(int j = 0;j < DV_Qv;++j)
			{
				P_InletShadowCA[k].fBP[i][j] = P_InletShadowCA[k].Cell_C[0]->fEq[i][j];
				P_InletShadowCA[k].fBP_x[i][j] = 0;
				P_InletShadowCA[k].fBP_y[i][j] = 0;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
				P_InletShadowCA[k].gBP[i][j] = P_InletShadowCA[k].Cell_C[0]->gEq[i][j];
				P_InletShadowCA[k].gBP_x[i][j] = 0;
				P_InletShadowCA[k].gBP_y[i][j] = 0;
#endif			
			}
		}
		#endif
		#ifdef _P_OUTLET_5_BCS_FLIP
			for(int k = 0;k < P_OutletFaceNum;++k)
			{
				for(int i = 0;i < DV_Qu;++i)
				for(int j = 0;j < DV_Qv;++j)
				{
					P_OutletShadowCA[k].fBP[i][j] = P_OutletShadowCA[k].Cell_C[0]->fEq[i][j];
					P_OutletShadowCA[k].fBP_x[i][j] = 0;
					P_OutletShadowCA[k].fBP_y[i][j] = 0;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
					P_OutletShadowCA[k].gBP[i][j] = P_OutletShadowCA[k].Cell_C[0]->gEq[i][j];
					P_OutletShadowCA[k].gBP_x[i][j] = 0;
					P_OutletShadowCA[k].gBP_y[i][j] = 0;
#endif					
				}
			}
		#endif
}
void ShockStructure()
{
	for(int n = 0;n < Cells;++n)
	{
		if(CellArray[n].xc <= 0)
		{
			CellArray[n].Rho = Rho0;
			CellArray[n].U = U0;
			CellArray[n].V = V0;
			CellArray[n].T = T0;
			CellArray[n].p = Rho0*R0*T0;
			CellArray[n].Lambda = 0.5/(T0*R0);
			CellArray[n].Mu = Mu0;
		}
		else
		{
			CellArray[n].Rho = Rho_Outlet;
			CellArray[n].U = U_Outlet;
			CellArray[n].V = V_Outlet;
			CellArray[n].T = T_Outlet;
			CellArray[n].p = Rho_Outlet*R0*T_Outlet;
			CellArray[n].Lambda = 0.5/(T_Outlet*R0);
			CellArray[n].Mu = Mu0*pow(T_Outlet/T0,Omega0);
		}
		CellArray[n].qx = 0;
		CellArray[n].qy = 0;
		CellArray[n].Factor();
		Update_phi_Eq(CellArray[n]);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			CellArray[n].fT[i][j] = CellArray[n].fEq[i][j];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
			CellArray[n].gT[i][j] = CellArray[n].gEq[i][j];
#endif
		}
	}
#ifdef _P_INLET_4_BCS_FLIP
	for(int k = 0;k < P_InletFaceNum;++k)
	{
		P_InletFaceA[k]->Rho_h = Rho0;
		P_InletFaceA[k]->U_h = U0;
		P_InletFaceA[k]->V_h = V0;
		P_InletFaceA[k]->T_h = T0;
		P_InletFaceA[k]->Lambda_h = 0.5/(T0*R0);
		P_InletFaceA[k]->p_h = Rho0*R0*T0;
		P_InletFaceA[k]->Mu_h = Mu0;
		P_InletFaceA[k]->qx_h = 0;
		P_InletFaceA[k]->qy_h = 0;
		P_InletFaceA[k]->Factor();
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_InletShadowCA[k].fBP[i][j] = P_InletFaceA[k]->lhsCell->fEq[i][j];
			P_InletShadowCA[k].fBP_x[i][j] = 0;		
			P_InletShadowCA[k].fBP_y[i][j] = 0;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
			P_InletShadowCA[k].gBP[i][j] = P_InletFaceA[k]->lhsCell->gEq[i][j];
			P_InletShadowCA[k].gBP_x[i][j] = 0;
			P_InletShadowCA[k].gBP_y[i][j] = 0;
#endif
		}
	}
#endif
#ifdef _P_OUTLET_5_BCS_FLIP
	for(int k = 0;k < P_OutletFaceNum;++k)
	{
		P_OutletFaceA[k]->Rho_h = Rho_Outlet;
		P_OutletFaceA[k]->U_h = U_Outlet;
		P_OutletFaceA[k]->V_h = V_Outlet;
		P_OutletFaceA[k]->T_h = T_Outlet;
		P_OutletFaceA[k]->Lambda_h = 0.5/(T_Outlet*R0);
		P_OutletFaceA[k]->p_h = Rho_Outlet*R0*T_Outlet;
		P_OutletFaceA[k]->Mu_h = Mu0*pow(T_Outlet/T0,Omega0);
		P_OutletFaceA[k]->qx_h = 0;
		P_OutletFaceA[k]->qy_h = 0;
		P_OutletFaceA[k]->Factor();
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_OutletShadowCA[k].fBP[i][j] = P_OutletFaceA[k]->lhsCell->fEq[i][j];
			P_OutletShadowCA[k].fBP_x[i][j] = 0;
			P_OutletShadowCA[k].fBP_y[i][j] = 0;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
			P_OutletShadowCA[k].gBP[i][j] = P_OutletFaceA[k]->lhsCell->gEq[i][j];
			P_OutletShadowCA[k].gBP_x[i][j] = 0;
			P_OutletShadowCA[k].gBP_y[i][j] = 0;
#endif
		}
	}
#endif
}
void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p)
{
	/*u = -U0*cos(2.0*PI*x)*sin(2.0*PI*y)*exp(-8.0*PI*PI*nu*t);
	v =  U0*sin(2.0*PI*x)*cos(2.0*PI*y)*exp(-8.0*PI*PI*nu*t);
	p = -0.25*U0*U0*(cos(4.0*PI*x) + cos(4.0*PI*y))*exp(-16.0*PI*PI*nu*t);*/
	u = -U0*cos(2.0*PI*x)*sin(2.0*PI*y)*exp(-8.0*PI*PI*Nu0*t)/(2.0*PI);
	v =  U0*sin(2.0*PI*x)*cos(2.0*PI*y)*exp(-8.0*PI*PI*Nu0*t)/(2.0*PI);
	p = -0.25*U0*U0*(cos(4.0*PI*x) + cos(4.0*PI*y))*exp(-16.0*PI*PI*Nu0*t)/(4.0*PI*PI);
}

/*void TG_Initialization()
{	
	for(int i = 0;i != Cells;++i)
	{
		double &x = CellArray[i].xc, &y = CellArray[i].yc;
		TaylorGreenVortex(0.0,CellArray[i].xc,CellArray[i].yc,
							  CellArray[i].u,CellArray[i].v,CellArray[i].p);
		CellArray[i].rho = Rho0 + CellArray[i].p/RT;
		Update_fEq(CellArray[i]);
//---------------------------------Initialize_fT-------------------------------
		double grad_ux = U0*sin(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vx = U0*cos(2.0*PI*x)*cos(2.0*PI*y);
//
		double grad_uy = -grad_vx;
		double grad_vy = -grad_ux;

		double grad_ut = nu*4.0*PI*U0*cos(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vt = -nu*4.0*PI*U0*sin(2.0*PI*x)*cos(2.0*PI*y);
		for(int k = 0;k < Q;++k)
		{
			double u1 = xi[k].u * CellArray[i].u + xi[k].v * CellArray[i].v;
			double A_u = (xi[k].u + u1*xi[k].u/RT - CellArray[i].u);
			double A_v = (xi[k].v + u1*xi[k].v/RT - CellArray[i].v);
			double A = omega[k]*Rho0/RT;
			double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
			double fEq_x = A * (A_u*grad_ux + A_v*grad_vx) * xi[k].u;
			double fEq_y = A * (A_u*grad_uy + A_v*grad_vy) * xi[k].v;
			double f = CellArray[i].fEq[k] - tau*(fEq_t + fEq_x + fEq_y);
			CellArray[i].fT[k] = f; //- 0.5*dt*(fEq_t + fEq_x + fEq_y);
		}
	}
}*/
void LidDrivenSquare()
{
	for(int n = 0;n < Cells;++n)
	{
		CellArray[n].Rho = Rho0;
		CellArray[n].U   = 0;
		CellArray[n].V   = 0;
		CellArray[n].T   = T0;
		CellArray[n].p   = Rho0*R0*T0;
		CellArray[n].Mu  = Mu0;
		CellArray[n].Lambda = Lambda0;
		CellArray[n].Factor();
		Update_phi_Eq(CellArray[n]);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			CellArray[n].fT[i][j] = CellArray[n].fEq[i][j];
			#ifndef _ARK_ISOTHERMAL_FLIP
			CellArray[n].gT[i][j] = CellArray[n].gEq[i][j];
			#endif
		}
	}
	for(int n = 0;n < Faces;++n)
	{
		FaceArray[n].Rho_h = Rho0;	
		if(fabs(Y_End - FaceArray[n].yf) < 1.0E-12)
		{
			FaceArray[n].U_h = U0;
			FaceArray[n].V_h = V0;
			FaceArray[n].rhsCell->U = U0;
			FaceArray[n].rhsCell->V = V0;
		}
		else
		{
			FaceArray[n].U_h = 0;
			FaceArray[n].V_h = 0;
			FaceArray[n].rhsCell->U = 0;
			FaceArray[n].rhsCell->V = 0;
		}
		FaceArray[n].T_h = T0;
		FaceArray[n].p_h = Rho0*R0*T0;
		FaceArray[n].Lambda_h = Lambda0;
		FaceArray[n].Mu_h = Mu0;
		FaceArray[n].Factor();
	}
	for(int k = 0;k < WallFaceNum;++k)
	{
		WallShadowCA[k].Rho = Rho0;
		WallShadowCA[k].T   = T0;
		WallShadowCA[k].Lambda = Lambda0;
		WallShadowCA[k].Mu = Mu0;
		WallShadowCA[k].Factor();
	}
}
/*void SquareInitialization()
{
	for(int i = 0;i != WallFaceNum;++i)
	{
		WallFaceA[i]->rhoh = Rho0;
		WallFaceA[i]->rhsCell->rho = Rho0;
		if(VelocityBCs == WallFaceA[i]->zone)
		{
			WallFaceA[i]->uh = U0;
			WallFaceA[i]->vh = 0.0;
			WallFaceA[i]->rhsCell->u = U0;
			WallFaceA[i]->rhsCell->v = 0.0;
			WallFaceA[i]->lhsCell->u = U0;
			WallFaceA[i]->lhsCell->v = 0.0;
		}
		else
		{
			WallFaceA[i]->uh = 0.0;
			WallFaceA[i]->vh = 0.0;
			WallFaceA[i]->rhsCell->u = 0.0;
			WallFaceA[i]->rhsCell->v = 0.0;
		}
	}
	for(int i = 0;i != Cells;++i)
	{
		CellArray[i].rho = Rho0;
		//CellArray[i].u = 0.0;
		//CellArray[i].v = 0.0;
		Update_phi_Eq(CellArray[i]);
		for(int k = 0;k < Q;++k)
			CellArray[i].fT[k] = CellArray[i].fEq[k];
	}
}*/
void TaylorCouetteAnalyticalSolution(double x,double y,double &u_A)
{
	const double eta = TC_r/TC_R;
	const double A = -W_i*eta*eta/(1.0 - eta*eta), B = W_i*TC_r*TC_r/(1.0 - eta*eta);
	double r = sqrt(x*x + y*y);
	u_A = A*r + B/r;	
}
/*void TaylorCouetteInitialization()
{
	for(int i = 0;i != WallFaceNum;++i)
	{
		Face_2D &face = *WallFaceA[i];
		WallFaceA[i]->rhoh = Rho0;
		WallFaceA[i]->rhsCell->rho = Rho0;
		if(VelocityBCs == WallFaceA[i]->zone)
		{
			WallFaceA[i]->uh = -U0*face.Vy;
			WallFaceA[i]->vh = U0*face.Vx;
			WallFaceA[i]->rhsCell->u = WallFaceA[i]->uh;
			WallFaceA[i]->rhsCell->v = WallFaceA[i]->vh;
		}
		else
		{
			WallFaceA[i]->uh = 0.0;
			WallFaceA[i]->vh = 0.0;
			WallFaceA[i]->rhsCell->u = 0.0;
			WallFaceA[i]->rhsCell->v = 0.0;
		}
	}
	for(int i = 0;i != Cells;++i)
	{
		CellArray[i].rho = Rho0;
		//CellArray[i].u = 0.0;
		//CellArray[i].v = 0.0;
		Update_fEq(CellArray[i]);
		for(int k = 0;k < Q;++k)
			CellArray[i].fT[k] = CellArray[i].fEq[k];
	}
}*/
void Riemann2D()
{
	for(int n = 0;n < Cells;++n)
	{
		if(CellArray[n].xc > 0.5 && CellArray[n].yc > 0.5)
		{
			CellArray[n].Rho = Rho1;
			CellArray[n].U = U1;
			CellArray[n].V = V1;
			CellArray[n].T = T1;
			CellArray[n].p = Rho1*R0*T1;
			CellArray[n].Lambda = 0.5/(T1*R0);
			CellArray[n].Mu = Mu0;
		}
		else if(CellArray[n].xc <= 0.5 && CellArray[n].yc > 0.5)
		{
			CellArray[n].Rho = Rho2;
			CellArray[n].U = U2;
			CellArray[n].V = V2;
			CellArray[n].T = T2;
			CellArray[n].p = Rho2*R0*T2;
			CellArray[n].Lambda = 0.5/(T2*R0);
			CellArray[n].Mu = Mu0*pow(T2/T0,Omega0);
		}
		else if(CellArray[n].xc <= 0.5 && CellArray[n].yc <= 0.5)
		{
			CellArray[n].Rho = Rho3;
			CellArray[n].U = U3;
			CellArray[n].V = V3;
			CellArray[n].T = T3;
			CellArray[n].p = Rho3*R0*T3;
			CellArray[n].Lambda = 0.5/(T3*R0);
			CellArray[n].Mu = Mu0*pow(T3/T0,Omega0);
		}
		else if(CellArray[n].xc > 0.5 && CellArray[n].yc <= 0.5)
		{
			CellArray[n].Rho = Rho4;
			CellArray[n].U = U4;
			CellArray[n].V = V4;
			CellArray[n].T = T4;
			CellArray[n].p = Rho4*R0*T4;
			CellArray[n].Lambda = 0.5/(T4*R0);
			CellArray[n].Mu = Mu0*pow(T4/T0,Omega0);
		}
		CellArray[n].qx = 0;
		CellArray[n].qy = 0;
		CellArray[n].Factor();
		Update_phi_Eq(CellArray[n]);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			CellArray[n].fT[i][j] = CellArray[n].fEq[i][j];
			CellArray[n].gT[i][j] = CellArray[n].gEq[i][j];
		}
	}
#ifdef _P_INLET_4_BCS_FLIP
	for(int k = 0;k < P_InletFaceNum;++k)
	{
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_InletShadowCA[k].fBP[i][j] = P_InletShadowCA[k].Cell_C[0]->fEq[i][j];
			P_InletShadowCA[k].gBP[i][j] = P_InletShadowCA[k].Cell_C[0]->gEq[i][j];
			P_InletShadowCA[k].fBP_x[i][j] = 0;
			P_InletShadowCA[k].gBP_x[i][j] = 0;
			P_InletShadowCA[k].fBP_y[i][j] = 0;
			P_InletShadowCA[k].gBP_y[i][j] = 0;
		}
	}
#endif
#ifdef _P_OUTLET_5_BCS_FLIP
	for(int k = 0;k < P_OutletFaceNum;++k)
	{
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_OutletShadowCA[k].fBP[i][j] = P_OutletShadowCA[k].Cell_C[0]->fEq[i][j];
			P_OutletShadowCA[k].gBP[i][j] = P_OutletShadowCA[k].Cell_C[0]->gEq[i][j];
			P_OutletShadowCA[k].fBP_x[i][j] = 0;
			P_OutletShadowCA[k].gBP_x[i][j] = 0;
			P_OutletShadowCA[k].fBP_y[i][j] = 0;
			P_OutletShadowCA[k].gBP_y[i][j] = 0;
		}
	}
#endif
}
