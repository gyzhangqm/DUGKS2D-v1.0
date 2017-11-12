#ifndef _MESH_2D_H_
#define _MESH_2D_H_

#include <cmath>
#include "../Header/ZeroConstant.h"

struct Cell_2D;

struct Face_2D;

struct Face_2D{
//
	Face_2D() = default;
	Face_2D(const Face_2D &rhs) = delete;
	Face_2D& operator=(const Face_2D &rhs) = delete;
	~Face_2D();
	Face_2D *shadowF = nullptr;
	int bc_type = 2;
	unsigned zone = 0;
	double *NodeX_F[2] = {nullptr,nullptr};
	double *NodeY_F[2] = {nullptr,nullptr};
	Cell_2D *lhsCell = nullptr,*rhsCell = nullptr;
	double xf = 0.0, yf = 0.0, Area = 0.0;
	double Vx = 0.0,Vy = 0.0;
	void SetArea();
	void SetNormalV();
//----------------------------------------------
	double **xi_n_dS = nullptr;
//
	double **fh    = nullptr;
	double **fBh   = nullptr;
	double **fEqh  = nullptr;
//
	double **gh    = nullptr;
	double **gBh   = nullptr;
	double **gEqh  = nullptr;
	static int const Qu = DV_Qu;
	static int const Qv = DV_Qv;
	void AllocateInFace();
//-----------macro variables---------------------
	double Rho_h = 0.0,U_h = 0.0,V_h = 0.0,p_h = 0.0,T_h = 0.0;
	double Lambda_h = 0.0,E_h = 0.0,qx_h = 0.0,qy_h = 0.0;
	double Mu_h = 0.0,Tau_h = 0.0;
//-----------------Factor------------------------
	double ah = 0.0, bh = 0.0, ch = 0.0,aQh = 0.0;
	void Factor();
private:
	int *use = new int(1);
};

struct Cell_2D{
//
	Cell_2D() = default;
	Cell_2D(const Cell_2D &rhs);
	Cell_2D& operator=(const Cell_2D &rhs);
	~Cell_2D();

//-----------Mesh------------------
	Cell_2D *ShadowC = nullptr;
	int celltype = 0;
	double xc = 0.0, yc = 0.0, volume = 0.0, DtSlashVolume = 0.0;
	double *NodeX_C[4] = {nullptr,nullptr,nullptr,nullptr};
	double *NodeY_C[4] = {nullptr,nullptr,nullptr,nullptr}; 
	Cell_2D *Cell_C[4] = {nullptr,nullptr,nullptr,nullptr};
	Face_2D *Face_C[4] = {nullptr,nullptr,nullptr,nullptr};
	int signFlux[4] = {0};
//---------------------------------------------------------
	double LS_M[2][2] = {{0.0,0.0},{0.0,0.0}};
	double wdx_C[4] = {0.0};
	double wdy_C[4] = {0.0};
	void SetVolume();
//----------probability density function---------------------
	double **fBP   = nullptr;				//	\bar{f}^+
	double **fBP_x = nullptr;			//	\frac{\partial {\bar{f}^+}}{\partial x}
	double **fBP_y = nullptr;			//  \frac{\partial {\bar{f}^+}}{\partial y}	
	double **fT    = nullptr;				//	\tilde{f}
	double **fEq   = nullptr;				//	f^{eq}
//
	double **gBP   = nullptr;
	double **gBP_x = nullptr;
	double **gBP_y = nullptr;
	double **gT    = nullptr;
	double **gEq   = nullptr;
	static int const Qu = DV_Qu;
	static int const Qv = DV_Qv;
	void AllocateInCell();
//----------macro variables---------------------
	double fFluxSum = 0, gFluxSum = 0;
	double U = 0.0,V = 0.0, Rho = 0.0, p = 0.0,E = 0.0,qx = 0.0,qy = 0.0;
	double shearTau[2][2] =  {{0.0,0.0},{0.0,0.0}};
	double Lambda = 0.0,T = 0.0,Tau = 0.0,Mu = 0.0;
	double U_1k = 0.0, V_1k = 0.0,T_1K = 0.0;
//-------------Factor---------------------
	double aBP = 0.0, bBP = 0.0, cBP = 0.0,aNEq = 0.0,aQ = 0.0;
	double fBPLimiter = 1,gBPLimiter = 1;
	void Factor();
private:
	int *use = new int(1);
};

inline 
void SetZero(double &d)
{
	if(d < infinitesimal && d > -infinitesimal) d= 0.0;
}
#endif