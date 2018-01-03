#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include "DUGKSDeclaration.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::ostringstream;
//--------------------------------------------------------
int Faces = 0,Nodes = 0,Cells = 0;

double *NodeX, *NodeY, *NodeZ;

Face_2D *FaceArray = nullptr;
Cell_2D *CellArray = nullptr;
//-------------------------------------------------------
//
//------------------------------------------------------
int InteriorFaceNum = 0, PeriodicFaceNum = 0, 
	WallFaceNum = 0, P_InletFaceNum = 0, P_OutletFaceNum = 0,
	SymmetryFaceNum = 0,P_FarfieldFaceNum = 0,V_InletFaceNum = 0, BoundFaceNum = 0;

Face_2D **InteriorFaceA = nullptr, **WallFaceA = nullptr,**PeriodicFaceA = nullptr,
		**P_InletFaceA = nullptr, **P_OutletFaceA = nullptr,**SymmetryFaceA = nullptr,
		**P_FarfieldFaceA = nullptr,**V_InletFaceA = nullptr,**BoundFaceA = nullptr;

Cell_2D *PeriodicShadowCA = nullptr, *WallShadowCA = nullptr,*P_InletShadowCA = nullptr,
		*P_OutletShadowCA = nullptr,*SymmetryShadowCA = nullptr,*P_FarfieldShadowCA = nullptr,
		*V_InletShadowCA = nullptr;

int step = 0;
double SumRho = 0.0, SumT = 0.0;
//--------------MeshConstruct.cpp----------------
extern int MeshConstruct(const string& s);
extern int MeshArea();
extern void FacesClassify();
extern int ShadowCellConstruct();
extern void NeighbourCellConstruct();
extern int MeshCheck();
extern int MeshOutput(const string& s);
extern void  Grad_LSMatrix();
//--------------Qmodel.cpp-------------------
extern void DiscreteVelocityAssign();
extern void setXiDotdS();
//--------------main.cpp---------------------
extern void AllocateResource();
extern void DeallocateResource();
//-------------Preprocess.cpp-------------------
extern void TG_Initialization();
extern void TaylorCouetteInitialization();
extern void ShockStructure();
extern void UniformFlow();
extern void Riemann2D();
extern void LidDrivenSquare();
//----------------Inc_2DSolver.cpp----------------
extern void DUGKS2DSolver();
//--------------Output.cpp------------------------
extern void OutputCase();
extern void Output_Convergence();
extern void SelfCheck();
int main()
{
	ostringstream oss_MeshName;
	oss_MeshName <<NL<< _MESHFILE_NAME_ARK;
	string MeshName(oss_MeshName.str().c_str());
	cout <<MeshName<<endl;
	OutputCase();
	SelfCheck();
	Output_Convergence();
//- - - - - - - - -Mesh- - - - - - - - - - - - -
	MeshConstruct(MeshName);
	MeshArea();
	FacesClassify();
	ShadowCellConstruct();
	NeighbourCellConstruct();
	#ifdef _ZERO_NDEBUG_FLIP
	MeshCheck();
	MeshOutput(MeshName);
	#endif
//-----------------Preprocess--------------------
	AllocateResource();
	Grad_LSMatrix();
//-------------------Initialization-------------------
	//UniformFlow();
	//LidDrivenSquare();
	//ShockStructure();
	//Riemann2D();
	//TG_Initialization();	
	//SquareInitialization();
	TaylorCouetteInitialization();
//------------------Solve-------------------
	#ifndef _ZERO_NDEBUG_FLIP
	DUGKS2DSolver();
	#endif
//------------------Afterprocess----------------
	//Output_uvProfile(TimeEnd);
	DeallocateResource();
	#ifdef _ZERO_NDEBUG_FLIP
	getchar();
	#endif
}
void SelfCheck()
{
//
	#ifdef _BB_BOUNDARY_SCHEME_FLIP
	if("BB" != _BC_ARK)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"\"BB\" != _BC_ARK"<<endl;
		getchar();
	}
	#endif
//
	#ifdef _NEE_BOUNDARY_SCHEME_FLIP
	if("NEE" != _BC_ARK)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"\"NEE\" != _BC_ARK"<<endl;
		getchar();
	}
	#endif
//
	#ifdef _CARTESIAN_MESH_FLIP
		if("CD" != _FLUX_SCHEME_ARK)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"CD\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	#else
		if("UW" != _FLUX_SCHEME_ARK)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"UW\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	#endif
	if("Quad" == _MESHTYPE_ARK || "Tri" == _MESHTYPE_ARK)
	{
		if("UW" != _FLUX_SCHEME_ARK)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"UW\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	}
	else if("Car" == _MESHTYPE_ARK)
	{
		if("CD" != _FLUX_SCHEME_ARK)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"CD\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	}
	else
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"Unknown Mesh Type"<<endl;
		getchar();
	}
}
//
void AllocateInCells(const int& CellNum,Cell_2D* &ptrCell)
{
	if(CellNum == 0) return;
	for(int k = 0;k < CellNum;++k)
		ptrCell[k].AllocateInCell();
}
void AllocateResource()
{
	DiscreteVelocityAssign();
//
	cout <<"Allocating Resources for Faces..."<<endl;
	for(int n = 0;n < Faces;++n)
		FaceArray[n].AllocateInFace();
	cout <<"Allocating Resources for Faces Done"<<endl;
//
	cout <<"Calculating xi*n*dS..."<<endl;
	setXiDotdS();
	cout <<"Calculating xi*n*dS Done"<<endl;
//
	cout <<"Allocating Resources for Cells..."<<endl;
	AllocateInCells(Cells,CellArray);
	AllocateInCells(WallFaceNum,WallShadowCA);
	AllocateInCells(PeriodicFaceNum,PeriodicShadowCA);
	AllocateInCells(P_InletFaceNum,P_InletShadowCA);
	AllocateInCells(P_OutletFaceNum,P_OutletShadowCA);
	AllocateInCells(SymmetryFaceNum,SymmetryShadowCA);
	AllocateInCells(P_FarfieldFaceNum,P_FarfieldShadowCA);
	AllocateInCells(V_InletFaceNum,V_InletShadowCA);
	for(int n = 0;n < PeriodicFaceNum;++n)
	{
		PeriodicShadowCA[n] = *PeriodicShadowCA[n].ShadowC;
	}
	cout <<"Allocating Resources for Cells Done"<<endl;
}
void DeallocateFaces(int BoundFaceNum, Face_2D** &ptrBoundFaceA)
{
	if (0 == BoundFaceNum) return;
	delete[] ptrBoundFaceA;
}
void DeallocateCells(int BoundFaceNum, Cell_2D* &ptrBoundShadowCA)
{
	if (0 == BoundFaceNum) return;
	delete[] ptrBoundShadowCA;
}
void DeallocateResource()
{
	delete[] xi_u;
	delete[] xi_v;
	delete[] NodeX;
	delete[] NodeY;
	delete[] FaceArray;
	delete[] CellArray;
	DeallocateFaces(InteriorFaceNum,InteriorFaceA);
	DeallocateFaces(PeriodicFaceNum,PeriodicFaceA);
	DeallocateFaces(WallFaceNum,WallFaceA);
	DeallocateFaces(P_InletFaceNum,P_InletFaceA);
	DeallocateFaces(P_OutletFaceNum,P_OutletFaceA);
	DeallocateFaces(SymmetryFaceNum,SymmetryFaceA);
	DeallocateFaces(P_FarfieldFaceNum,P_FarfieldFaceA);
	DeallocateFaces(V_InletFaceNum,V_InletFaceA);
//
	DeallocateCells(WallFaceNum,WallShadowCA);
	DeallocateCells(P_InletFaceNum,P_InletShadowCA);
	DeallocateCells(P_OutletFaceNum,P_OutletShadowCA);
	DeallocateCells(SymmetryFaceNum,SymmetryShadowCA);
	DeallocateCells(P_FarfieldFaceNum,P_FarfieldShadowCA);
	DeallocateCells(V_InletFaceNum,V_InletShadowCA);
	//DeallocateCells(PeriodicFaceNum,PeriodicShadowCA);
}