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
//--------------Preprocess.cpp-------------------
extern void AllocateResource();
extern void DeallocateResource();
extern void TG_Initialization();
extern void SquareInitialization();
extern void TaylorCouetteInitialization();
extern void ShockStructure();
extern void UniformFlow();
extern void Riemann2D();
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
	ShockStructure();
	//Riemann2D();
	//TG_Initialization();	
	//SquareInitialization();
	//TaylorCouetteInitialization();
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
/*int main()
{
	clock_t c_start,c_end;
	time_t t_start,t_end;
	c_start = clock();
	time(&t_start);

	c_end = clock();
	time(&t_end);
	//cout << "time used : " << difftime(c_end,c_start)/CLOCKS_PER_SEC << "s";
	//cout << "time used : " << t_end - t_start << "s";

}*/

