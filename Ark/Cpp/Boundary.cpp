#include "DUGKSDeclaration.h"
#include <iostream>
using std::cout;

void P_Inlet_4_Boundary()
{
	for(int k = 0;k < P_InletFaceNum;++k)
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_InletShadowCA[k].fBP[i][j] = P_InletShadowCA[k].Cell_C[0]->fBP[i][j];
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			P_InletShadowCA[k].gBP[i][j] = P_InletShadowCA[k].Cell_C[0]->gBP[i][j];
			#endif
		}
}
void P_Outlet_5_Boundary()
{
	for(int k = 0;k < P_OutletFaceNum;++k)
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_OutletShadowCA[k].fBP[i][j] = P_OutletShadowCA[k].Cell_C[0]->fBP[i][j];
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			P_OutletShadowCA[k].gBP[i][j] = P_OutletShadowCA[k].Cell_C[0]->gBP[i][j];
			#endif
		}
}
void WallShadowC_fBP(Cell_2D &shadowCell)
{
// used for GradfBP
		Cell_2D const *cell = shadowCell.Cell_C[0];
		shadowCell.Rho = cell->Rho;//Non-Equilibrium Extrapolation
		Update_phi_Eq(shadowCell);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			shadowCell.fBP[i][j] = shadowCell.fEq[i][j]
			+ cell->aBP*(cell->fT[i][j] - cell->fEq[i][j]);
	//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP	
			shadowCell.gBP[i][j] = shadowCell.gEq[i][j]
			+ cell->aBP*(cell->gT[i][j] - cell->gEq[i][j]);
		#endif
		}
}
// void Wall_3_Boundary(Face_2D &face)
// {
// 	face.Rho_h = face.lhsCell->Rho;
// 	Update_phi_Eqh(face);
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.fh[i][j] = face.fEqh[i][j]
// 		+ face.lhsCell->aNEq*(face.lhsCell->fT[i][j] - face.lhsCell->fEq[i][j]);
// 		face.fh[i][j] *= face.xi_n_dS[i][j];
// 	//isothermal flip
// 	#ifndef _ARK_ISOTHERMAL_FLIP	
// 		face.gh[i][j] = face.gEqh[i][j]
// 		+ face.lhsCell->aNEq*(face.lhsCell->gT[i][j] - face.lhsCell->gEq[i][j]);
// 		face.gh[i][j] *= face.xi_n_dS[i][j];
// 	#endif
// 	}
// }
