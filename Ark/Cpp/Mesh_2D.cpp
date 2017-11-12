#include "Mesh_2D.h"

double TriArea(double Xbeg, double Ybeg, double X_A, double Y_A,double X_B, double Y_B)
{
	double Area = 0.5 * ((X_A - Xbeg) * (Y_B - Ybeg) - (X_B - Xbeg) * (Y_A - Ybeg));
	return (Area > 0 ? Area : -Area);
}
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

Cell_2D::Cell_2D(const Cell_2D &rhs)
{
	fBP   = rhs.fBP;
	fBP_x = rhs.fBP_x;
	fBP_y = rhs.fBP_y;
	fT    = rhs.fT;
	fEq   = rhs.fEq;
//
	gBP   = rhs.gBP;
	gBP_x = rhs.gBP_x;
	gBP_y = rhs.gBP_y;
	gT    = rhs.gT;
	gEq   = rhs.gEq;
	use   = rhs.use;
//
	++*rhs.use;
}
Cell_2D& Cell_2D::operator=(const Cell_2D &rhs)
{
	++*rhs.use;
	if(--*use == 0)
	{	
		DeallocateARK(fBP,Qu,Qv);
		DeallocateARK(fBP_x,Qu,Qv);
		DeallocateARK(fBP_y,Qu,Qv);
		DeallocateARK(fT,Qu,Qv);
		DeallocateARK(fEq,Qu,Qv);
	//
#ifndef _ARK_ISOTHERMAL_FLIP
		DeallocateARK(gBP,Qu,Qv);
		DeallocateARK(gBP_x,Qu,Qv);
		DeallocateARK(gBP_y,Qu,Qv);
		DeallocateARK(gT,Qu,Qv);
		DeallocateARK(gEq,Qu,Qv);
#endif
		delete use;
	}
	fBP   = rhs.fBP;
	fBP_x = rhs.fBP_x;
	fBP_y = rhs.fBP_y;
	fT    = rhs.fT;
	fEq   = rhs.fEq;
//
	gBP   = rhs.gBP;
	gBP_x = rhs.gBP_x;
	gBP_y = rhs.gBP_y;
	gT    = rhs.gT;
	gEq   = rhs.gEq;
//
	use   = rhs.use;
	return *this;
}
Cell_2D::~Cell_2D()
{
	if(--*use == 0)
	{	
		DeallocateARK(fBP,Qu,Qv);
		DeallocateARK(fBP_x,Qu,Qv);
		DeallocateARK(fBP_y,Qu,Qv);
		DeallocateARK(fT,Qu,Qv);
		DeallocateARK(fEq,Qu,Qv);
//
#ifndef _ARK_ISOTHERMAL_FLIP
		DeallocateARK(gBP,Qu,Qv);
		DeallocateARK(gBP_x,Qu,Qv);
		DeallocateARK(gBP_y,Qu,Qv);
		DeallocateARK(gT,Qu,Qv);
		DeallocateARK(gEq,Qu,Qv);
#endif
//
		delete use;
	}
}
void Cell_2D::AllocateInCell()
{
	AllocateARK(fBP,Qu,Qv);
	AllocateARK(fBP_x,Qu,Qv);
	AllocateARK(fBP_y,Qu,Qv);
	AllocateARK(fT,Qu,Qv);
	AllocateARK(fEq,Qu,Qv);
//
#ifndef _ARK_ISOTHERMAL_FLIP
	AllocateARK(gBP,Qu,Qv);
	AllocateARK(gBP_x,Qu,Qv);
	AllocateARK(gBP_y,Qu,Qv);
	AllocateARK(gT,Qu,Qv);
	AllocateARK(gEq,Qu,Qv);
#endif
//
	use = new int(1);
}
void Cell_2D::Factor()
{
	Tau = Mu/p;
	aBP = (2.0*Tau - h)/(2.0*Tau + dt);
	bBP = 1.0 - aBP;
	cBP = Tau*bBP;
	aNEq = 2.0*Tau/(2.0*Tau + dt);
	aQ = Tau/(2.0*Tau + dt*Pr);
}
void Cell_2D::SetVolume()
{
	if(3 == celltype) 
	{
		xc = (*NodeX_C[0] + *NodeX_C[1] + *NodeX_C[2])/3.0;
		yc = (*NodeY_C[0] + *NodeY_C[1] + *NodeY_C[2])/3.0;
		volume = TriArea(*NodeX_C[0],*NodeY_C[0],*NodeX_C[1],*NodeY_C[1],*NodeX_C[2],*NodeY_C[2]);
	}
	else if(4 == celltype)
	{
//		
		double xc_A,yc_A,xc_B,yc_B,volume_A,volume_B;
//		
		xc_A = (*NodeX_C[0] + *NodeX_C[1] + *NodeX_C[2])/3.0;
		yc_A = (*NodeY_C[0] + *NodeY_C[1] + *NodeY_C[2])/3.0;
		volume_A = TriArea(*NodeX_C[0],*NodeY_C[0],*NodeX_C[1],*NodeY_C[1],*NodeX_C[2],*NodeY_C[2]);
//
		xc_B = (*NodeX_C[0] + *NodeX_C[3] + *NodeX_C[2])/3.0;
		yc_B = (*NodeY_C[0] + *NodeY_C[3] + *NodeY_C[2])/3.0;
		volume_B = TriArea(*NodeX_C[0],*NodeY_C[0],*NodeX_C[3],*NodeY_C[3],*NodeX_C[2],*NodeY_C[2]);
//		
		volume = volume_A + volume_B;
//
		xc = (volume_A*xc_A + volume_B*xc_B)/volume;
		yc = (volume_A*yc_A + volume_B*yc_B)/volume;
	}
	DtSlashVolume = dt/volume;
}
//
Face_2D::~Face_2D()
{
//
	if(--*use == 0)	
	{
		DeallocateARK(fh,Qu,Qv);
		DeallocateARK(fBh,Qu,Qv);
		DeallocateARK(fEqh,Qu,Qv);
//
#ifndef _ARK_ISOTHERMAL_FLIP
		DeallocateARK(gh,Qu,Qv);
		DeallocateARK(gBh,Qu,Qv);
		DeallocateARK(gEqh,Qu,Qv);
#endif
		DeallocateARK(xi_n_dS,Qu,Qv);
//
		delete use;
	}
}
void Face_2D::Factor()
{
	Tau_h = Mu_h/p_h;
	ah = 2.0*Tau_h/(2.0*Tau_h + h);
	bh = 1.0 - ah;
	ch = Tau_h*bh;
	aQh = Tau_h/(2.0*Tau_h + h*Pr);
}
void Face_2D::AllocateInFace()
{
	AllocateARK(fh,Qu,Qv);
	AllocateARK(fBh,Qu,Qv);
	AllocateARK(fEqh,Qu,Qv);
//
#ifndef _ARK_ISOTHERMAL_FLIP
	AllocateARK(gh,Qu,Qv);
	AllocateARK(gBh,Qu,Qv);
	AllocateARK(gEqh,Qu,Qv);
#endif
	AllocateARK(xi_n_dS,Qu,Qv);
//
	use = new int(1);
}
void Face_2D::SetArea()
{
	xf = 0.5*(*(NodeX_F[0]) + *(NodeX_F[1]));
	yf = 0.5*(*(NodeY_F[0]) + *(NodeY_F[1]));
	Area = sqrt(
				(*NodeY_F[1] - *NodeY_F[0]) * (*NodeY_F[1] - *NodeY_F[0]) +
				(*NodeX_F[1] - *NodeX_F[0]) * (*NodeX_F[1] - *NodeX_F[0])
			   );
}
void Face_2D::SetNormalV()
{
	double dy = (*NodeY_F[1] - *NodeY_F[0]), 
	 	   dx = (*NodeX_F[1] - *NodeX_F[0]);
	//SetZero(dx);SetZero(dy); 	 	   
	Vx = dy/Area;Vy = -dx/Area;
}