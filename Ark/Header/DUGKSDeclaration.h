#ifndef _DUGKS_DECLARATION_H_
#define _DUGKS_DECLARATION_H_

#include "Mesh_2D.h"

//--------------------------------------------------
extern int Faces,Nodes,Cells;
extern double *NodeX, *NodeY, *NodeZ;
extern Face_2D *FaceArray;
extern Cell_2D *CellArray;
//--------------------------------------------------
//
//--------------------------------------------------
extern int InteriorFaceNum,PeriodicFaceNum,WallFaceNum,P_InletFaceNum,
		   P_OutletFaceNum,	SymmetryFaceNum, P_FarfieldFaceNum, V_InletFaceNum,BoundFaceNum;

extern Face_2D **InteriorFaceA,**WallFaceA,**PeriodicFaceA,**P_InletFaceA,
	           **P_OutletFaceA,**SymmetryFaceA,**P_FarfieldFaceA,**V_InletFaceA,**BoundFaceA;

extern Cell_2D *PeriodicShadowCA, *WallShadowCA,*P_InletShadowCA,*P_OutletShadowCA,
			   *SymmetryShadowCA, *P_FarfieldShadowCA,*V_InletShadowCA;

extern Cell_2D *PeriodicShadowC_NE, *PeriodicShadowC_NW,
			   *PeriodicShadowC_SE, *PeriodicShadowC_SW;

extern double * const xi_u, * const xi_v;
//---------------------------------------------------

extern int step;

extern double SumRho, SumT;
//------------------------function declaration---------------------

//----------DmVn.cpp----------------------------
extern void Update_phi_Eq(Cell_2D &cell);

extern void Update_phi_Eqh(Face_2D &face);

extern void Update_MacroVar(Cell_2D& cell);

extern void Update_MacroVar_h(Face_2D& face);

// extern void CD_Interior_phi_Bh(Face_2D &face,int i,int j);

// extern void UW_Interior_phi_Bh(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j);

// extern void UW_Interior_phi_Bh_Limiter(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j);

// extern void VenkatakrishnanFluxLimiter(Cell_2D &cell,int const &i,int const &j);

extern void Update_phi_Eqh(Face_2D &face,int i,int j);

extern void Update_phi_h(Face_2D& face,int i,int j);

extern void Update_phiFlux_h(Face_2D &face,int i,int j);

//
template<typename T>
inline int MeshIndex(const T &End,const T &Beg)
{
	return (nullptr == End ? 0 : End - Beg + 1);
}

#endif
