#include <iostream>
#include <iomanip>
#include <omp.h>
#include "DUGKSDeclaration.h"
using std::ios;
using std::setiosflags;
using std::setprecision;
using std::cout;
using std::endl;

double const aTP = 4.0/3.0, bTP = 1.0/3.0;

int const ThreadNum = omp_get_max_threads();

double ResidualPer1k = 1.0;

double const PRECISION = 1.0E-16;

double const XDEBUG = 13.25,YDEBUG = 0.25;


//--------------------------DEBUG-------------------------------

extern void Output_L2Norm(double const &t,double &L2_uv, double &L2_p);

extern void Output_Flowfield(double const &t,int step);

extern void Output_SumRho(double t);

extern void Output_Residual(double t,double Residual);

extern void Output_UVP(double const &t);

//------------------------Boundary.cpp----------------------

extern void P_Inlet_4_Boundary();

extern void P_Outlet_5_Boundary();

extern void WallShadowC_fBP(Cell_2D &shadowCell);

extern void Wall_3_Boundary(Face_2D &face);

//----------------------------------DEBUG---------------------------------------
extern void Output_fBh(Face_2D& face,double t);

extern void Output_gBh(Face_2D& face,double t);

extern void Output_fh(Face_2D& face,double t);

extern void Output_gh(Face_2D& face,double t);

extern void Output_fT(Cell_2D& face,double t);

extern void Output_gT(Cell_2D& face,double t);

extern void Output_xcyc();

extern void Output_fT_Append(Cell_2D &cell,double dt);

extern void Output_gT_Append(Cell_2D &cell,double dt);

extern void Output_fh_Append(Face_2D &face,double dt);

extern void Output_gh_Append(Face_2D &face,double dt);

// extern void Output_fFlux_Append(Cell_2D &cell,double dt);

// extern void Output_gFlux_Append(Cell_2D &cell,double dt);

extern void Output_phi_Bh(Face_2D &face,double t);

//------------------------------------------------------------------------------

void LeastSquareDebug();

void LeastSquareABC(Cell_2D *center);

void Update_phi_BP(Cell_2D& cell);

void Update_phi_h(Face_2D& face);

void Update_phiFlux_h(Face_2D& face);

void Update_Flux(Face_2D &face);

void Update_BoundFlux(Face_2D &face);

void Update_phi_T(Cell_2D& cell);

void Zero_PartialDerivatives(Cell_2D& cell);

void Flux_2D(Face_2D &face);

void Flux_2D_Limiter(Face_2D &face);

void Update_Residual(int step);

void UpdateL2Error(int step);

void Update_SumRho(int step);
//--------------------------------DmVn.cpp--------------------------
extern void IntegralShearStress();

void DUGKS2DSolver()
{
Output_Flowfield(0.0,0);
cout << "ThreadNum : "<<ThreadNum<<endl;
omp_set_num_threads(ThreadNum);
step = 0;
#pragma omp parallel
{
//while(step < 10000)
while(ResidualPer1k > RESIDUAL)
{
	#pragma omp for schedule(guided)
	for(int n = 0;n < Cells;++n)
	{
		Update_phi_Eq(CellArray[n]);
		Update_phi_BP(CellArray[n]);
	}
//--------------------------------------Update-fBP--------------------------------
	#ifdef _Wall_3_BCs_FLIP
	#pragma omp for schedule(guided)
	for(int n = 0;n < WallFaceNum;++n)
		WallShadowC_fBP(WallShadowCA[n]);
	#endif
	#ifdef _P_INLET_4_BCS_FLIP	
		P_Inlet_4_Boundary();
	#endif
	#ifdef _P_OUTLET_5_BCS_FLIP
		P_Outlet_5_Boundary();
	#endif
//-------------------------------LeastSquare------------------------------
	 #pragma omp for schedule(guided)
	 for(int n = 0;n < Cells;++n)
	 {
		LeastSquareABC(CellArray + n);
		//Zero_PartialDerivatives(CellArray[n]);
	 }
//-------------------------------Flux-------------------------------------
//-------------------------------Interior Face-----------------------------	
	#ifdef _ARK_LIMITER_FLIP
		#pragma omp for schedule(guided)
		for(int n = 0;n < InteriorFaceNum;++n)
		{
			Flux_2D_Limiter(*InteriorFaceA[n]);
		}
		#pragma omp for schedule(guided)
		for(int n = 0;n < BoundFaceNum;++n)
			Flux_2D(*BoundFaceA[n]);
	#else
		#pragma omp for schedule(guided)
		for(int n = 0;n < InteriorFaceNum;++n)
		{
			Flux_2D(*InteriorFaceA[n]);
		}
	#endif
	#ifdef _Wall_3_BCs_FLIP
		#pragma omp for schedule(guided)
		for(int n = 0;n < WallFaceNum;++n)
		{
			Wall_3_Boundary(*WallFaceA[n]);
		}
	#endif
// //----------------------------------------Wall Face----------------------------------
// 	#ifdef _Wall_3_BCs_FLIP
// // 	for(int i_Wall = 0;i_Wall < WallFaceNum;++i_Wall)
// // 	{
		
// // 		#ifdef _NEE_BOUNDARY_SCHEME_FLIP
// // 		NonEquilibriumExtrapolation(*WallFaceA[i_Wall]);
// // 		#endif
// // //		
// // 		//DiffusiveScatter(*WallFaceA[i_Wall]);
// // //		
// // 		#ifdef _BB_BOUNDARY_SCHEME_FLIP
// // 		BounceBack(*WallFaceA[i_Wall]);
// // 		#endif
// // 	}
// 	for(int i = 0;i < WallFaceNum;++i)
// 		Update_BoundFlux(*WallFaceA[i]);
// 	#endif
//----------------------------------------------------------------------------------
	#pragma omp for schedule(guided)
	for(int n = 0;n < Cells;++n)
		Update_phi_T(CellArray[n]);
//
	#pragma omp for schedule(guided)
	for(int n = 0;n < Cells;++n)
		 Update_MacroVar(CellArray[n]);
	#pragma omp single
	{
		++step;
		if(step%ConvergenceControl == 0)
		{
			Update_SumRho(step);
			#ifdef _OUTPUT_L2NORM_ERROR_FLIP
			UpdateL2Error(step);
			#endif
		}
		if(step%ResidualControl == 0)
			Update_Residual(step);
	}
}
}
	IntegralShearStress();
	Output_Flowfield((step)*dt,step);
}
void LeastSquareDebug()
{
	for(int i = 0;i < Cells;++i)
	{
		Cell_2D *center = CellArray + i, *neighbour = nullptr;
		for(int m = 0;m < DV_Qu;++m)
		for(int n = 0;n < DV_Qv;++n)
		{
			center->fBP[m][n] = -1.0;
			for(int Iface = 0;Iface < center->celltype;++Iface)
			{
			neighbour = center->Cell_C[Iface];
				if(neighbour->xc == center->xc && neighbour->yc > center->yc)
				{
					neighbour->fBP[m][n] = 1.5;
				}
				else if(neighbour->xc == center->xc && neighbour->yc < center->yc)
				{
					neighbour->fBP[m][n] = 0.0;
				}
				else if(neighbour->yc == center->yc && neighbour->xc > center->xc)
				{
					neighbour->fBP[m][n] = 5;
				}
				else if(neighbour->yc == center->yc && neighbour->xc < center->xc)
				{
					neighbour->fBP[m][n] = 2;
				}
				else
				{
					cout << "neither xc nor yc is the same;"<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
				}
			}
			double Sum_wdxdfBP = 0.0;
			double Sum_wdydfBP = 0.0;
			for(int Iface = 0;Iface < center->celltype;++Iface)
			{
				neighbour = center->Cell_C[Iface];
				Sum_wdxdfBP += center->wdx_C[Iface]*(neighbour->fBP[m][n] - center->fBP[m][n]);
				Sum_wdydfBP += center->wdy_C[Iface]*(neighbour->fBP[m][n] - center->fBP[m][n]);
			}
			center->fBP_x[m][n] = center->LS_M[0][0]*Sum_wdxdfBP + center->LS_M[0][1]*Sum_wdydfBP;
			center->fBP_y[m][n] = center->LS_M[1][0]*Sum_wdxdfBP + center->LS_M[1][1]*Sum_wdydfBP;
			if(CellArray[i].fBP_x[m][n] != 3.0*0.5*NL || CellArray[i].fBP_y[m][n] != 1.5*0.5*NL)
			{
				cout <<"CellArray : "<<i<<" ---------"<<endl;
				cout <<"m = "<<m <<"  "<<CellArray[i].fBP_x[m][n]
					<<"    "<<CellArray[i].fBP_y[m][n]<<endl;
				getchar();
			}
		}
	}
}
void LeastSquareABC(Cell_2D *center)
{
	Cell_2D  *neighbour = nullptr;
	for(int m = 0;m < DV_Qu;++m)
	for(int n = 0;n < DV_Qv;++n)
	{
		double Sum_wdxdfBP = 0.0;
		double Sum_wdydfBP = 0.0;
		double Sum_wdxdgBP = 0.0;
		double Sum_wdydgBP = 0.0;
		for(int Iface = 0;Iface < center->celltype;++Iface)
		{
			neighbour = center->Cell_C[Iface];
			Sum_wdxdfBP += center->wdx_C[Iface]*(neighbour->fBP[m][n] - center->fBP[m][n]);
			Sum_wdydfBP += center->wdy_C[Iface]*(neighbour->fBP[m][n] - center->fBP[m][n]);
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			Sum_wdxdgBP += center->wdx_C[Iface]*(neighbour->gBP[m][n] - center->gBP[m][n]);
			Sum_wdydgBP += center->wdy_C[Iface]*(neighbour->gBP[m][n] - center->gBP[m][n]);
			#endif
		}
		center->fBP_x[m][n] = center->LS_M[0][0]*Sum_wdxdfBP + center->LS_M[0][1]*Sum_wdydfBP;
		center->fBP_y[m][n] = center->LS_M[1][0]*Sum_wdxdfBP + center->LS_M[1][1]*Sum_wdydfBP;
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		center->gBP_x[m][n] = center->LS_M[0][0]*Sum_wdxdgBP + center->LS_M[0][1]*Sum_wdydgBP;
		center->gBP_y[m][n] = center->LS_M[1][0]*Sum_wdxdgBP + center->LS_M[1][1]*Sum_wdydgBP;
		#endif
	}
}
void Zero_PartialDerivatives(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		cell.fBP_x[i][j] = 0.0;
		cell.fBP_y[i][j] = 0.0;
//isothermal flip
#ifndef _ARK_ISOTHERMAL_FLIP
		cell.gBP_x[i][j] = 0.0;
		cell.gBP_y[i][j] = 0.0;
#endif
	}
}
void Update_phi_BP(Cell_2D& cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		cell.fBP[i][j] = cell.aBP*cell.fT[i][j] + cell.bBP*cell.fEq[i][j];
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		cell.gBP[i][j] = cell.aBP*cell.gT[i][j] + cell.bBP*cell.gEq[i][j];
		#endif
	}
}
void Update_phi_h(Face_2D& face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		face.fh[i][j] = face.ah*face.fBh[i][j] + face.bh*face.fEqh[i][j];
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		face.gh[i][j] = face.ah*face.gBh[i][j] + face.bh*face.gEqh[i][j];
		#endif
	}
}
void Update_phiFlux_h(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		face.fh[i][j] *= face.xi_n_dS[i][j];
//isothermal flip
#ifndef _ARK_ISOTHERMAL_FLIP
		face.gh[i][j] *= face.xi_n_dS[i][j];
#endif
	}
}
inline
double VenkatakrishnanExpression(double a,double b)
{
	return (a*a + 2*a*b + Kh3)/(a*a + 2*b*b + a*b + Kh3);
} 
void VenkatakrishnanFluxLimiter(Cell_2D &cell,int const &i,int const &j)
{
	double GradfBPDotDelta[4],LfBP[4];
	double GradgBPDotDelta[4],LgBP[4];
//
	double MaxfBP = -1E+5,MinfBP = -MaxfBP, MaxgBP = -1E+5,MingBP = -MaxgBP;
	double MinfBPLimiter = 1,MingBPLimiter = 1;
//	
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		if(cell.Cell_C[iFace]->fBP[i][j] > MaxfBP)
			MaxfBP = cell.Cell_C[iFace]->fBP[i][j];
		if(cell.Cell_C[iFace]->fBP[i][j] < MinfBP)
			MinfBP = cell.Cell_C[iFace]->fBP[i][j];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
		if(cell.Cell_C[iFace]->gBP[i][j] > MaxgBP)
			MaxgBP = cell.Cell_C[iFace]->gBP[i][j];
		if(cell.Cell_C[iFace]->gBP[i][j] < MingBP)
			MingBP = cell.Cell_C[iFace]->gBP[i][j];
#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		GradfBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.fBP_x[i][j]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.fBP_y[i][j];					        
		if(GradfBPDotDelta[iFace] > 0)
			LfBP[iFace] =  VenkatakrishnanExpression(MaxfBP-cell.fBP[i][j],GradfBPDotDelta[iFace]);
		else if(GradfBPDotDelta[iFace] < 0)
			LfBP[iFace] =  VenkatakrishnanExpression(MinfBP-cell.fBP[i][j],GradfBPDotDelta[iFace]);
		else
			LfBP[iFace] = 1;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
		GradgBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.gBP_x[i][j]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.gBP_y[i][j];
		if(GradgBPDotDelta[iFace] > 0)
			LgBP[iFace] =  VenkatakrishnanExpression(MaxgBP-cell.gBP[i][j],GradgBPDotDelta[iFace]);
		else if(GradgBPDotDelta[iFace] < 0)
			LgBP[iFace] =  VenkatakrishnanExpression(MingBP-cell.gBP[i][j],GradgBPDotDelta[iFace]);
		else
			LgBP[iFace] = 1;
#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		if(LfBP[iFace] < MinfBPLimiter) MinfBPLimiter = LfBP[iFace];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
		if(LgBP[iFace] < MingBPLimiter) MingBPLimiter = LgBP[iFace];
#endif
	}
	cell.fBPLimiter = MinfBPLimiter;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
	cell.gBPLimiter = MingBPLimiter;
#endif
}
/*void DiffusiveScatter(Face_2D &face)
{
	double uu = face.uh*face.uh + face.vh*face.vh;
	double u1,xi_dot_n;
	double rhohWall_A = 0.0, rhohWall_b = 0.0;
	for(int k = 0;k < Q;++k)
	{
		xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		u1 = xi[k].u*face.uh + xi[k].v*face.vh;
		if(xi_dot_n > 0)
		{
			Solve_Interior_fBh(face,face.lhsCell,k);
			rhohWall_b += xi_dot_n*(face.ah*face.fBh[k] + face.bh*ConstanceInfEq(u1,uu,k));
			rhohWall_A += xi_dot_n*bh*omega[k];
		}
		else if(xi_dot_n < 0)
		{
			rhohWall_b += xi_dot_n*ConstanceInfEq(u1,uu,k);
			rhohWall_A += xi_dot_n*omega[k];
		}
	}
	face.rhoh = -rhohWall_b/rhohWall_A;
	Update_fEqh(face);
	double out = 0.0,in = 0.0;
	for(int k = 0;k < Q;++k)
	{
		xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		if(xi_dot_n > 0)
		{
			face.fh[k] = face.ah*face.fBh[k] + face.bh*face.fEqh[k];
			out += face.fh[k]*xi_dot_n;
		}
		else if(xi_dot_n < 0)
		{
			face.fh[k] = face.fEqh[k];
			in += face.fh[k]*xi_dot_n;
		}
	}
	if(in + out > 1.0E-16 || in + out < -1.0E-16)
	{
		cout <<in + out<<endl;
		getchar();
	}
}

void NonEquilibriumExtrapolation(Face_2D &face)
{
	Cell_2D &cell = *face.lhsCell;
//
	face.rhoh = cell.rho;
//
	Update_fEqh(face);
	for(int k = 0;k != Q;++k)
	{
		face.fh[k] =  face.fEqh[k] + cell.aNEq*(cell.fT[k] - cell.fEq[k]);
	}
}
void BounceBack(Face_2D &face)
{
	face.rhoh = face.lhsCell->rho;
	Update_fEqh(face);
	for(int k = 0;k < Q;++k)
	{
		double xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		if(xi_dot_n > 0.0)
		{
			Solve_Interior_fBh(face,face.lhsCell,k);
		}
	}
	for(int k = 0;k < Q;++k)
	{
		double xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		if(xi_dot_n < 0.0)
		{
			double xi_dot_u = face.uh * xi[k].u + face.vh * xi[k].v;
			if(k < 5)
			{	
				face.fBh[k] = face.fBh[k + 4] + 2.0*face.rhoh*omega[k]*xi_dot_u/RT;
			}
			else
			{
				face.fBh[k] = face.fBh[k - 4] + 2.0*face.rhoh*omega[k]*xi_dot_u/RT;
			}

		}
	}
	Update_fh(face);
}*/
void Update_phi_fBh(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#ifdef _CARTESIAN_MESH_FLIP
			CD_Interior_phi_Bh(face,i,j);
		#else
		if(face.xi_n_dS[i][j] >= 0)
			UW_Interior_phi_Bh(face,face.lhsCell,i,j);
		else
			UW_Interior_phi_Bh(face,face.rhsCell,i,j);
		#endif
	}
}
void Update_phi_fBh_Limiter(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#ifdef _CARTESIAN_MESH_FLIP
			CD_Interior_phi_Bh(face,i,j);
		#else
		if(face.xi_n_dS[i][j] >= 0)
			UW_Interior_phi_Bh_Limiter(face,face.lhsCell,i,j);
		else
			UW_Interior_phi_Bh_Limiter(face,face.rhsCell,i,j);
		#endif
	}
}
void Flux_2D_Limiter(Face_2D &face)
{	
	Update_phi_fBh_Limiter(face);
	Update_MacroVar_h(face);
	Update_phi_Eqh(face);
	Update_phi_h(face);
	Update_phiFlux_h(face);
//------------------------------------------DEBUG---------------------------------
}
void Flux_2D(Face_2D &face)
{	
	Update_phi_fBh(face);
	Update_MacroVar_h(face);
	Update_phi_Eqh(face);
	Update_phi_h(face);
	Update_phiFlux_h(face);
}
//-------------------------------------------------------------------------------
void Update_phi_T(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		//cell.fT[i][j] = aTP*cell.fBP[i][j] - bTP*cell.fT[i][j] + cell.DtSlashVolume *cell.fFlux[i][j];
		//cell.gT[i][j] = aTP*cell.gBP[i][j] - bTP*cell.gT[i][j] + cell.DtSlashVolume *cell.gFlux[i][j];
		cell.fFluxSum = 0;
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		cell.gFluxSum = 0;
		#endif
		for(int k = 0;k < cell.celltype;++k)
		{
			cell.fFluxSum += cell.signFlux[k]*cell.Face_C[k]->fh[i][j];
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			cell.gFluxSum += cell.signFlux[k]*cell.Face_C[k]->gh[i][j];
			#endif
		}
		cell.fT[i][j] = aTP*cell.fBP[i][j] - bTP*cell.fT[i][j] + cell.DtSlashVolume *cell.fFluxSum;
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		cell.gT[i][j] = aTP*cell.gBP[i][j] - bTP*cell.gT[i][j] + cell.DtSlashVolume *cell.gFluxSum;
		#endif
	}
}
void Update_SumRho(int step)
{
	SumRho = 0.0;
	SumT = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		SumRho += CellArray[i].Rho;
		#ifndef _ARK_ISOTHERMAL_FLIP 
		SumT += CellArray[i].T;
		#endif
	}
	Output_SumRho(step*dt);
}
void Update_Residual(int step)
{
		double SumUV = 0.0,Sumdudv = 0.0;
		double du = 0.0, dv = 0.0;
		for(int i = 0;i < Cells;++i)
		{
			du = CellArray[i].U - CellArray[i].U_1k;
			dv = CellArray[i].V - CellArray[i].V_1k;
			Sumdudv += du*du + dv*dv;
			SumUV += CellArray[i].U*CellArray[i].U + CellArray[i].V*CellArray[i].V;
			CellArray[i].U_1k = CellArray[i].U;
			CellArray[i].V_1k = CellArray[i].V;
		}
		ResidualPer1k = sqrt(Sumdudv/(SumUV + 1.0E-30));
		Output_Residual(step*dt,ResidualPer1k);
		if(step%writeFileControl == 0)
		Output_Flowfield(step*dt, step);
	#ifndef _ARK_NOHUP_FLIP
		cout << setiosflags(ios::scientific) << setprecision(6);
		cout <<step <<"    "<<step*dt<<"    "<<SumRho<<"    "<<SumT<<"    "<<ResidualPer1k<<'\n';
	#endif
}
void UpdateL2Error(int step)
{
	double L2_uv, L2_p;
	Output_L2Norm(step*dt,L2_uv,L2_p);
	//cout << setiosflags(ios::scientific) << setprecision(12);
	//cout <<step <<"    "<<SumRho<<"    "<<L2_uv<<"    "<<L2_p<<'\n';
}