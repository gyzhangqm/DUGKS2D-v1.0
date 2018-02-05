#include "DUGKSDeclaration.h"

namespace D2Q9{

double const xi_u[DV_Qv] = {0,1,1,0,-1,-1,-1,0,1};

double const xi_v[DV_Qv] = {0,0,1,1,1,0,-1,-1,-1};

}

extern double * const xi_u = new double[DV_Qv];

extern double * const xi_v = new double[DV_Qv];

const double omega[DV_Qv]={4.0/9.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0};

void DiscreteVelocityAssign()
{
	for(int j = 0;j < DV_Qv;++j)
	{
		xi_u[j] = MaSpan*D2Q9::xi_u[j];
		xi_v[j] = MaSpan*D2Q9::xi_v[j];
	}
}
// void setXiDotdS()
// {
// 	for(int n = 0;n < Faces;++n)
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		FaceArray[n].xi_n_dS[i][j] = 
// 		FaceArray[n].Area*(xi_u[j]*FaceArray[n].Vx + xi_v[j]*FaceArray[n].Vy);
// 	}
// }
void Update_phi_Eq(Cell_2D &cell)
{
	double uu,u1,u2;
	uu = cell.U*cell.U + cell.V*cell.V;
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		u1 = (xi_u[j]*cell.U + xi_v[j]*cell.V)/RT;
		cell.fEq[i][j] = omega[j]*(cell.Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
	}
}
void Update_phi_Eqh(Face_2D &face)
{
	double uu,u1,u2;
	uu = face.U_h*face.U_h + face.V_h*face.V_h;
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		u1 = (xi_u[j]*face.U_h + xi_v[j]*face.V_h)/RT; 
		face.fEqh[i][j] = omega[j]*(face.Rho_h + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
	}
}
double vectorDot(const double *const first,const double *const second)
{
	double sum = 0;
	for(int j = 0;j < DV_Qv;++j)
		sum += first[j]*second[j];
	return sum;
}
double vectorDot(const double *const first)
{
	double sum = 0;
	for(int j = 0;j < DV_Qv;++j)
		sum += first[j];
	return sum;
}
void Update_MacroVar(Cell_2D& cell)
{
	cell.Rho  = vectorDot(cell.fT[0]);
	cell.U    = (vectorDot(cell.fT[0],xi_u) + h*cell.Fx)/Rho0;
	cell.V    = (vectorDot(cell.fT[0],xi_v) + h*cell.Fy)/Rho0;
	//cell.U    = (vectorDot(cell.fT[0],xi_u))/Rho0;
	//cell.V    = (vectorDot(cell.fT[0],xi_v))/Rho0;
	cell.p    = 0.5*(cell.Rho - Rho0)/Lambda0;
}
void Update_MacroVar_h(Face_2D& face)
{
	face.Rho_h  = vectorDot(face.fBh[0]);
	face.U_h    = (vectorDot(face.fBh[0],xi_u) + 0.5*h*face.Fx_h)/Rho0;
	face.V_h    = (vectorDot(face.fBh[0],xi_v) + 0.5*h*face.Fy_h)/Rho0;
	//face.U_h    = (vectorDot(face.fBh[0],xi_u))/Rho0;
	//face.V_h    = (vectorDot(face.fBh[0],xi_v))/Rho0;
	// face.p_h    = 0.5*(face.Rho_h - Rho0)/Lambda0;
}
// void UW_Interior_phi_Bh(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
// {
// 	double dx = face.xf - h*xi_u[j] - ptr_C->xc;
// 	double dy = face.yf - h*xi_v[j] - ptr_C->yc;
// 	face.fBh[i][j] = ptr_C->fBP[i][j] + (dx*ptr_C->fBP_x[i][j] + dy*ptr_C->fBP_y[i][j]);
// }
// void CD_Interior_phi_Bh(Face_2D &face,int i,int j)
// {
// 	double _dx = face.lhsCell->xc - face.rhsCell->xc;
// 	double _dy = face.lhsCell->yc - face.rhsCell->yc;
// 	SetZero(_dx);
// 	SetZero(_dy);
// 	double _fBP_xF,_fBP_yF;
// 	{
// 		if(0.0 == _dx)
// 		{
// 			 _fBP_xF =  0.5*(face.lhsCell->fBP_x[i][j] + face.rhsCell->fBP_x[i][j]);
// 			 _fBP_yF = (face.lhsCell->fBP[i][j] - face.rhsCell->fBP[i][j])/_dy;
// 		}
// 		else if(0.0 == _dy)
// 		{
// 			_fBP_yF = 0.5*(face.lhsCell->fBP_y[i][j] + face.rhsCell->fBP_y[i][j]);
// 			_fBP_xF = (face.lhsCell->fBP[i][j] - face.rhsCell->fBP[i][j])/_dx;
// 		}
// 		face.fBh[i][j] = 0.5*(face.lhsCell->fBP[i][j] + face.rhsCell->fBP[i][j])
// 				 - h*(_fBP_xF*xi_u[j] + _fBP_yF*xi_v[j]);
// 	}	
// }
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	cell.force[i][j] = (cell.Fx*(xi_u[j]-cell.U) + cell.Fy*(xi_v[j]-cell.V))
// 						*2*Lambda0*cell.fEq[i][j]/Rho0;
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	face.forceh[i][j] = (face.Fx_h*(xi_u[j]-face.U_h)+face.Fy_h*(xi_v[j]-face.V_h))
// 						*2*Lambda0*face.fEqh[i][j]/Rho0;
// 	}
// }
// void Wall_3_Boundary(Face_2D &face)
// {

// }
void IntegralShearStress(){}
//
//-----------------------------Inc unsteady Taylor Green vortex------------------
extern void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p);
void unsteadyTaylorGreen()
{	
	for(int n = 0;n != Cells;++n)
	{
		double x = CellArray[n].xc, y = CellArray[n].yc;
		TaylorGreenVortex(0.0,x,y,CellArray[n].U,CellArray[n].V,CellArray[n].p);
		CellArray[n].Rho = Rho0 + CellArray[n].p/RT;
		CellArray[n].T = T0;
		CellArray[n].Lambda = Lambda0;
		CellArray[n].Mu = Mu0;
		CellArray[n].Tau = 2.0*Nu0*Lambda0;
		CellArray[n].Factor();
		Update_phi_Eq(CellArray[n]);
//---------------------------------Initialize_fT-------------------------------
		double grad_ux = U0*sin(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vx = U0*cos(2.0*PI*x)*cos(2.0*PI*y);
//
		double grad_uy = -grad_vx;
		double grad_vy = -grad_ux;

		double grad_ut = Nu0*4.0*PI*U0*cos(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vt = -Nu0*4.0*PI*U0*sin(2.0*PI*x)*cos(2.0*PI*y);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			double u1 = (xi_u[j]) * CellArray[n].U + (xi_v[j]) * CellArray[n].V;
			double A_u = ((xi_u[j]) + u1*(xi_u[j])/RT - CellArray[n].U);
			double A_v = ((xi_v[j]) + u1*(xi_v[j])/RT - CellArray[n].V);
			double A = omega[j]*Rho0/RT;
			double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
			double fEq_x = A * (A_u*grad_ux + A_v*grad_vx) * (xi_u[j]);
			double fEq_y = A * (A_u*grad_uy + A_v*grad_vy) * (xi_v[j]);
			double f = CellArray[n].fEq[i][j] - CellArray[n].Tau*(fEq_t + fEq_x + fEq_y);
			CellArray[n].fT[i][j] = f; //- 0.5*dt*(fEq_t + fEq_x + fEq_y);
		}
	}
	for(int n = 0;n < Faces;++n)
	{
		FaceArray[n].T_h = T0;
		FaceArray[n].Mu_h = Mu0;
		FaceArray[n].Lambda_h = Lambda0;
		FaceArray[n].Tau_h = 2.0*Nu0*Lambda0;
		FaceArray[n].Factor();
	}
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