#ifndef _MESH_CONSTRUCTFUN_H_
#define _MESH_CONSTRUCTFUN_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "DUGKSDeclaration.h"

using std::string;
using std::ios;
using std::ofstream;
using std::istringstream;
using std::ostringstream;
using std::cout;
using std::endl;

const int MeshPerLine = 5,NumPerCell = 4;

void Allocate_Mesh(const int &index,const int &Num)
{
	if(index == 10)
	{		
		Nodes = Num;
		NodeX = new double[Nodes];
		NodeY = new double[Nodes];
	}
	else if(index == 12)
	{
		Cells = Num;
		CellArray = new Cell_2D[Cells]();
	}
	else if(index == 13)
	{
		Faces = Num;
		FaceArray = new Face_2D[Faces]();
	}
	else
	{
		cout << "unknown head type  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<<endl;
	}
}

void ConstructNodes(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	InFile_Mesh >> std::dec;
	for(int i = 0;i != Nodes;++i)
	{
		InFile_Mesh >> NodeX[i] >> NodeY[i];
		#ifdef _ZERO_NDEBUG_FLIP
		if(i%ZeroDebugControl == 0)
		std::cout <<std::setiosflags(ios::scientific)<<std::setprecision(15)
				  << NodeX[i] << "  " << NodeY[i] <<std::setprecision(6)
				  << resetiosflags(ios::scientific) << std::endl;
		#endif
	}
	cout <<"Nodes Constructed" << endl;
}
void ConstructCells(int* const &ptrHexLine,const int &count)
{
	if(5 != count) 
	{
		std::cerr <<"Error : " <<  __FILE__ <<" : in function "<< __func__
				  <<" at line " <<__LINE__ << endl;
		getchar();
		return;
	}
	if(3 == ptrHexLine[count - 1] )
	{
		for(int i = 0;i != Cells;++i)
			CellArray[i].celltype = 4;
			
	}
	else if(1 == ptrHexLine[count - 1])
	{
		for(int i = 0;i != Cells;++i)
			CellArray[i].celltype = 3;
	}
	else
	{
		std::cerr << "Error : " <<  __FILE__ <<" : in function "<< __func__
				  <<" at line " <<__LINE__ << endl;
		std::cerr << "wrong element type" << endl;
		getchar();
	}
}
void ConstructCells(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	int tmp;
	for(int i = 0;i != Cells;++i)
	{
		InFile_Mesh >> tmp;
		if(tmp == 1)
			CellArray[i].celltype = 3;
		else if(tmp == 3)
			CellArray[i].celltype = 4;
		else
		{
			cout << "wrong element index for 2D domain" << endl;
			getchar();
		}
	}
}
void PushNode(Cell_2D &Cell,const int &BegN,const int &EndN)
{
	if(nullptr == Cell.NodeX_C[0])
	{
		Cell.NodeX_C[0] = NodeX+BegN;
		Cell.NodeX_C[1] = NodeX+EndN;
		Cell.NodeY_C[0] = NodeY+BegN;
		Cell.NodeY_C[1] = NodeY+EndN;
	}
	else if(nullptr == Cell.NodeX_C[2] || nullptr == Cell.NodeX_C[3])
	{
		if((NodeX + BegN) == Cell.NodeX_C[1]  && (NodeY + BegN) == Cell.NodeY_C[1])
		{
			Cell.NodeX_C[2] = NodeX + EndN;
			Cell.NodeY_C[2] = NodeY + EndN;
		}
		else if((NodeX + EndN) == Cell.NodeX_C[1]  && (NodeY + EndN) == Cell.NodeY_C[1])
		{
			Cell.NodeX_C[2] = NodeX + BegN;
			Cell.NodeY_C[2] = NodeY + BegN;
		}
		else if((NodeX + BegN) == Cell.NodeX_C[0]  && (NodeY + BegN) == Cell.NodeY_C[0])
		{
			Cell.NodeX_C[3] = NodeX + EndN;
			Cell.NodeY_C[3] = NodeY + EndN;
		}
		else if((NodeX + EndN) == Cell.NodeX_C[0]  && (NodeY + EndN) == Cell.NodeY_C[0])
		{
			Cell.NodeX_C[3] = NodeX + BegN;
			Cell.NodeY_C[3] = NodeY + BegN;
		}
	}
}
void PushFace(Cell_2D &Cell,Face_2D *ptrF)
{
	int i;
	for(i = 0;Cell.Face_C[i] != nullptr;++i);
	Cell.Face_C[i] = ptrF;
}
void PushCell(Cell_2D &Cell,Cell_2D *ptrF)
{
	int i;
	for(i = 0;Cell.Cell_C[i] != nullptr;++i);
	Cell.Cell_C[i] = ptrF;
}
void ConstructFacesCells(const int &FC,const int (&N_C)[MeshPerLine],const int &bc_type,const int &zone);
void ConstructFaces(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	const int &zone = ptrHexLine[0];
	const int &Beg_Face = ptrHexLine[1] -1;
	const int &End_Face = ptrHexLine[2];
	const int &bc_type = ptrHexLine[3];
	const int &face_type = ptrHexLine[4];
	int N_C[MeshPerLine] = {0,0,0,0,0};//set of node index and cell index 

	if(0 == face_type)
	{		
		InFile_Mesh >> std::hex;
		for(int i = Beg_Face;i < End_Face;++i)
		{
			InFile_Mesh >> N_C[0] >> N_C[1] >> N_C[2] >> N_C[3] >> N_C[4];
			ConstructFacesCells(i,N_C,bc_type,zone);
		}
		InFile_Mesh >> std::dec;			
	}
	else if(2 == face_type)
	{
		InFile_Mesh >> std::hex;
		for(int i = Beg_Face;i < End_Face;++i)
		{
			InFile_Mesh >> N_C[1] >> N_C[2] >> N_C[3] >> N_C[4];
			ConstructFacesCells(i,N_C,bc_type,zone);
		}
		InFile_Mesh >> std::dec;
	}
	else
	{
		cout <<__FILE__ <<"  "<<__LINE__<<"  "<<__func__<<endl;
		getchar();
	}
}
void ConstructFacesCells(const int &FC,const int (&N_C)[MeshPerLine],const int &bc_type,const int &zone)
{
	int BegN = N_C[1], EndN = N_C[2], lhsC = N_C[3], rhsC = N_C[4];
	--BegN; --EndN; --lhsC; --rhsC;
	FaceArray[FC].bc_type = bc_type;
	FaceArray[FC].zone = zone;
//
	FaceArray[FC].NodeX_F[0] = NodeX+BegN;
	FaceArray[FC].NodeX_F[1] = NodeX+EndN;
//
	FaceArray[FC].NodeY_F[0] = NodeY+BegN;
	FaceArray[FC].NodeY_F[1] = NodeY+EndN;
	FaceArray[FC].lhsCell = CellArray + lhsC;
//	
	PushNode(CellArray[lhsC],BegN,EndN);
	PushFace(CellArray[lhsC],FaceArray+FC);

	if(bc_type == 2)
	{
		FaceArray[FC].rhsCell = CellArray + rhsC;
		PushNode(CellArray[rhsC],EndN,BegN);
		PushFace(CellArray[rhsC],FaceArray+FC);
		//PushCell(CellArray[rhsC],CellArray + lhsC);
		//PushCell(CellArray[lhsC],CellArray + rhsC);
		++InteriorFaceNum;
	}
	else if(12 == bc_type || 8 == bc_type)
	{
		++PeriodicFaceNum;
		++BoundFaceNum;
	}
	else if(3 == bc_type)
	{
		++WallFaceNum;
		++BoundFaceNum;
	}
	else if(4 == bc_type)
	{
		++P_InletFaceNum;
		++BoundFaceNum;
	}
	else if(5 == bc_type)
	{
		++P_OutletFaceNum;
		++BoundFaceNum;
	}
	else if(7 == bc_type)
	{
		++SymmetryFaceNum;
		++BoundFaceNum;
	}
	else if(9 == bc_type)
	{
		++P_FarfieldFaceNum;
		++BoundFaceNum;
	}
	else if(10 == bc_type)
	{
		++V_InletFaceNum;
		++BoundFaceNum;
	}
}
int MeshArea()
{
	double dL_Min = 1.0;
	cout << "Traversing Faces : " << endl;
	for(int i = 0;i < Faces;++i)
	{
		FaceArray[i].SetArea();
		FaceArray[i].SetNormalV();
		if(FaceArray[i].Area < dL_Min) dL_Min = FaceArray[i].Area;
		#ifdef _ZERO_NDEBUG_FLIP
		if(i%ZeroDebugControl == 0)
		cout <<"Face : "<<i<<" "<<FaceArray[i].bc_type<<"  "
			 <<std::setiosflags(ios::scientific)<<std::setprecision(6)
			 <<FaceArray[i].xf <<" "<<FaceArray[i].yf<<" "
		     <<FaceArray[i].Area <<" "<<FaceArray[i].Vx<<" "<<FaceArray[i].Vy
		     <<resetiosflags(ios::scientific)<<endl;
		#endif
	}
	if(fabs(dL_Min - MinL) > 1.0E-15)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<std::setiosflags(ios::scientific)<<std::setprecision(15)
		<<"Mesh MinL : "<<dL_Min<<"  Const MinL : "<<MinL<<"    "<<fabs(dL_Min - MinL)<<endl;
		ofstream OutFile_Case("../FlowField/global/dMinL.ark",ios::app);
		if(!OutFile_Case)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"dMinL.ark open failed."<<endl;
			getchar();
			return 0;
		}
		OutFile_Case <<std::setiosflags(ios::scientific)<<std::setprecision(15)
					 <<dL_Min<<"    =    dL_Min"<<endl;
		getchar();
	}
	cout << "SetArea Done, SetNormalV Done" <<endl;
	cout << "Traversing Cells : " << endl;
	for(int i = 0;i != Cells;++i)
	{
		CellArray[i].SetVolume();
		#ifdef _ZERO_NDEBUG_FLIP
		if(i%ZeroDebugControl == 0)
		cout <<"Cell : "<<i<<" "<<std::setiosflags(ios::scientific)<<std::setprecision(6)
		 	 <<CellArray[i].xc <<" "<<CellArray[i].yc<<" "<< CellArray[i].volume
		 	 <<resetiosflags(ios::scientific)<<endl;
		#endif
	}
	cout << "SetVolume Done" <<endl;
	return 0;
}
void PeriodFaces(const int &F0,const int &F1);
void ConstructPeriodicFaces(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	const int &Beg = ptrHexLine[0] -1;
	const int &End = ptrHexLine[1];
	int F0 = 0;
	int F1 = 0;
	InFile_Mesh >> std::hex;
	for(int i = Beg;i < End;++i)
	{
		InFile_Mesh >> F0 >> F1;
		--F0;--F1;
		PeriodFaces(F0,F1);
	}
	InFile_Mesh >> std::dec;
}
void PeriodFaces(const int &F0,const int &F1)
{
	if(nullptr == FaceArray[F0].rhsCell && nullptr == FaceArray[F1].rhsCell)
	{
		FaceArray[F0].shadowF = FaceArray + F1; 
		FaceArray[F1].shadowF = FaceArray + F0; 
	}
	else
	{
		cout << FaceArray[F0].rhsCell <<"  " << FaceArray[F1].rhsCell <<endl;
		cout <<__FILE__ <<"  "<<__LINE__<<"  "<<__func__<<endl;
		getchar();
	}
}

#endif
