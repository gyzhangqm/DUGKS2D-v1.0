#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include "MeshConstructFunc.h"

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

void EraseParenthese(string &string_line,int &body)
{
	while(' ' == string_line.back() || '\n' == string_line.back() || '\r' == string_line.back())
		string_line.pop_back();
	if(string_line.back() == '(') ++body;
	auto It = remove_if(string_line.begin(),string_line.end(),
									[](const char& c){return (c == '(' || c == ')');});
	if(It == string_line.end()) body = 99;
	string_line.erase(It,string_line.end());
}
int StringToHex(int* const &ptrHexLine, istringstream &iss_line,int &count)
{
	string str[MeshPerLine];
	for(;iss_line >> str[count];++count);
	if(0 == count) return 0;
	#ifdef _ZERO_NDEBUG_FLIP
		cout << count << endl;
		for(int i = 0;i != count;++i)
			cout << str[i] << " ";
		cout << endl;
	#endif

	for(int i = 0;i != count; ++i)
		ptrHexLine[i] = stoi(str[i],nullptr,16);
	#ifdef _ZERO_NDEBUG_FLIP
		for(int i = 0;i != count;++i)
			cout << ptrHexLine[i] << " ";
		cout << endl;
	#endif
	return 0;
}
void HeadProcess(const int &index,int* const &ptrHexLine)
{
	Allocate_Mesh(index,ptrHexLine[2]);
}
void BodyProcess(const int &index,ifstream& InFile_Mesh,int* const &ptrHexLine)
{
	if(index == 10)
	{		
		ConstructNodes(ptrHexLine,InFile_Mesh);
	}
	else if(index == 12)
	{
		ConstructCells(ptrHexLine,InFile_Mesh);
	}
	else if(index == 13)
	{
		ConstructFaces(ptrHexLine,InFile_Mesh);
	}
	else if(index == 18)
	{
		ConstructPeriodicFaces(ptrHexLine,InFile_Mesh);
	}
	else
	{
		cout << "Invalid index during body processing" << endl;
		getchar();
	}
}
int MeshConstruct(const string &s)
{
	int *ptrHexLine = new int[MeshPerLine];
	int line = 0, count = 0, index = 0,body = 0,FaceCount = 0;
	string string_line;
	ifstream InFile_Mesh;
	InFile_Mesh.open("/home/yangzr/Mesh/"+ s +".cas");
	if(!InFile_Mesh)
	{
		cout << "file open failed  " <<__FILE__ <<" : "<<__LINE__<<"  "<<__func__<<endl; 
		return -99999;
	}
	while(getline(InFile_Mesh,string_line))
	{
		if(string_line.size() == 0) continue;
		cout << string_line << endl;
		body = 0;index = 0;	count = 0;
		for(int i = 0;i != MeshPerLine;++i)
			ptrHexLine[i] =  0;
		EraseParenthese(string_line,body);
		if(99 == body) continue;		
		istringstream iss_line(string_line);
		if(!(iss_line >> index) || !(index == 10 || index == 12 || index == 13 || index == 18))
			continue;
		StringToHex(ptrHexLine,iss_line,count);
		if(ptrHexLine[0] == 0)
			HeadProcess(index,ptrHexLine);
		else
		{
			if(body)
			{
				BodyProcess(index,InFile_Mesh,ptrHexLine);
			}
			else if(12 == index)//
			{
				ConstructCells(ptrHexLine,count);
			}
		}
	}
	cout << "Faces " << Faces << " Nodes " << Nodes << " Cells " << Cells <<endl;
	delete []ptrHexLine;
	return 0;
}

int MeshOutput(const string& s)
{
	cout << "Mesh Output verifing..."<<endl;
	ofstream OutFile_Mesh;
	OutFile_Mesh.open("../MeshOutPut/"+ s +".plt");
	if(!OutFile_Mesh)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		return -99999;
	}
	ostringstream fname,zonename,dataNE;
	zonename<<"ZONE T = \"Mesh\"\n";
	dataNE<<"Nodes="<<Nodes<<", Elements="<<Cells<<", ZONETYPE=FEQuadrilateral\n";
	string tecformat[5]={"VARIABLES = \"X\",\"Y\"\n",
						zonename.str().c_str(),dataNE.str().c_str(),"DATAPACKING=BLOCK\n",
						"VarLocation=([1-2]=NODAL)\n"};
	OutFile_Mesh << tecformat[0]<<tecformat[1]<<tecformat[2]<<tecformat[3]<<tecformat[4];
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_Mesh << NodeX[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_Mesh << endl;
	}
	OutFile_Mesh << endl;
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_Mesh << NodeY[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_Mesh << endl;
	}
	OutFile_Mesh << endl;
	for(int i = 0;i != Cells;++i)
	{
		OutFile_Mesh << MeshIndex(CellArray[i].NodeX_C[0] , NodeX)<< " "
					 << MeshIndex(CellArray[i].NodeX_C[1] , NodeX)<< " "
					 << MeshIndex(CellArray[i].NodeX_C[2] , NodeX)<< " " 
					 << MeshIndex(CellArray[i].NodeX_C[3] , NodeX)<< endl;
	}
	cout <<"Mesh Output Done" <<endl;
	OutFile_Mesh.close();
	cout <<"Face Output verifing..."<<endl;
	ofstream OutFile_Face;
	OutFile_Face.open("../MeshOutPut/13_" + s + ".dat");
	if(!OutFile_Face)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
	}
	OutFile_Face << std::hex;
	for(int i = 0;i != InteriorFaceNum;++i)
	{
		OutFile_Face << MeshIndex(FaceArray[i].NodeY_F[0] , NodeY)<<" "
					 << MeshIndex(FaceArray[i].NodeY_F[1] , NodeY)<<" "
					 << MeshIndex(FaceArray[i].lhsCell , CellArray)<<" "
					 << MeshIndex(FaceArray[i].rhsCell , CellArray)<<endl;
	}
	OutFile_Face.close();
//
	OutFile_Face.open("../MeshOutPut/18_" + s +".dat");
	if(!OutFile_Face)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
	}
	OutFile_Face << std::hex;
	for(int i = 0;i != Faces;++i)
	{
		if(8 == FaceArray[i].bc_type)
		OutFile_Face << MeshIndex(FaceArray[i].shadowF , FaceArray)<<" "
					 << i + 1 <<endl;
	}
	OutFile_Face.close();
	cout <<"Face Output Done" <<endl;
	cout <<"Node Output verifing..."<<endl;
	ofstream OutFile_Node;
	OutFile_Node.open("../MeshOutPut/10_" + s + ".dat");
	if(!OutFile_Node)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
	}
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_Node <<std::setiosflags(ios::scientific) <<std::setprecision(16)
					 <<NodeX[i]<<"  "<<NodeY[i] << endl;
	}
	OutFile_Node.close();
	cout <<"Node Output Done"<<endl;
	return 0;
}
void AllocateFaces(int BoundFaceNum, Face_2D** &ptrBoundFace, vector<Face_2D*> &vecBoundFace)
{
	if (0 == BoundFaceNum) return;
	ptrBoundFace  = new Face_2D*[BoundFaceNum];
	for(int i = 0;i != BoundFaceNum;++i)
		ptrBoundFace[i] = vecBoundFace[i];
}
void FacesClassify()
{
	vector<Face_2D*> InteriorVec,WallVec,PeriodicVec,P_InletVec,P_OutletVec,V_InletVec,
					 SymmetryVec,P_FarfieldVec,BoundVec;//2,3,12 -8,4,5,10,7,9;
	for(int i = 0;i != Faces; ++i)
	{
		if(2 == FaceArray[i].bc_type)
			InteriorVec.push_back(&FaceArray[i]);
		else if(12 == FaceArray[i].bc_type || 8 == FaceArray[i].bc_type)
		{
			PeriodicVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(3 == FaceArray[i].bc_type)
		{
			WallVec.push_back(&FaceArray[i]);
			BoundVec.push_back(FaceArray+i);
		}
		else if(4 == FaceArray[i].bc_type)
		{
			P_InletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(5 == FaceArray[i].bc_type)
		{
			P_OutletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(7 == FaceArray[i].bc_type)
		{
			SymmetryVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(9 == FaceArray[i].bc_type)
		{
			P_FarfieldVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(10 == FaceArray[i].bc_type)
		{
			V_InletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else
		{
			cout << "FaceArray : " << i <<" unknown bc_type : " << FaceArray[i].bc_type
			<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
			getchar();
		}
	}
	#ifdef _ZERO_NDEBUG_FLIP
	if(InteriorFaceNum != InteriorVec.size())
	{
		cout << "FaceArray : " << InteriorFaceNum <<" interior Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(WallFaceNum != WallVec.size())
	{
		cout << "FaceArray : " << WallFaceNum <<" Wall Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(PeriodicFaceNum != PeriodicVec.size())
	{
		cout << "FaceArray : " << PeriodicFaceNum <<" Periodic Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_InletFaceNum != P_InletVec.size())
	{
		cout << "FaceArray : " << P_InletFaceNum <<" Pressure Inlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_OutletFaceNum != P_OutletVec.size())
	{
		cout << "FaceArray : " << P_OutletFaceNum <<" Pressure Outlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(V_InletFaceNum != V_InletVec.size())
	{
		cout << "FaceArray : " << V_InletFaceNum <<" Velocity Inlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(SymmetryVec.size() != SymmetryFaceNum)
	{
		cout << "FaceArray : " << SymmetryFaceNum <<" Symmetry Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_FarfieldVec.size() != P_FarfieldFaceNum)
	{
		cout << "FaceArray : " << P_FarfieldFaceNum <<" Pressure Farfield Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(BoundVec.size() != BoundFaceNum)
	{
		cout << "FaceArray : " << BoundFaceNum <<" Boundary Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	cout <<"Faces : "<<Faces<<'\n'
		 <<"BoundFaceNum : "<<BoundFaceNum<<'\n'
		 <<"InFaceNum : "<<InteriorFaceNum<<'\n'
		 <<"WallFaceNum : "<<WallFaceNum<<'\n'
		 <<"PeriodicFaceNum : "<<PeriodicFaceNum<<'\n'
		 <<"P_InletFaceNum : "<<P_InletFaceNum<<'\n'
		 <<"P_OutletFaceNum : "<<P_OutletFaceNum<<'\n'
		 <<"SymmetryFaceNum : "<<SymmetryFaceNum<<'\n'
		 <<"P_FarfieldFaceNum : "<<P_FarfieldFaceNum<<'\n'
		 <<"V_InletFaceNum : "<<V_InletFaceNum<<endl;
	if(InteriorFaceNum + WallFaceNum + PeriodicFaceNum + P_InletFaceNum + P_OutletFaceNum
		+ SymmetryFaceNum + P_FarfieldFaceNum + V_InletFaceNum != Faces)
	{
		_PRINT_ERROR_MSG_FLIP
		cout << "Faces != InteriorFaceNum + WallFaceNum"<<endl;
		getchar();
	}
	if(WallFaceNum + PeriodicFaceNum + P_InletFaceNum + P_OutletFaceNum
		+ SymmetryFaceNum + P_FarfieldFaceNum + V_InletFaceNum != BoundFaceNum)
	{
		_PRINT_ERROR_MSG_FLIP
		cout << "BoundFaceNum != ..."<<endl;
		getchar();
	}
	#endif
	AllocateFaces(InteriorFaceNum,InteriorFaceA,InteriorVec);
	AllocateFaces(PeriodicFaceNum,PeriodicFaceA,PeriodicVec);
	AllocateFaces(WallFaceNum,WallFaceA,WallVec);
	AllocateFaces(P_InletFaceNum,P_InletFaceA,P_InletVec);
	AllocateFaces(P_OutletFaceNum,P_OutletFaceA,P_OutletVec);
	AllocateFaces(SymmetryFaceNum,SymmetryFaceA,SymmetryVec);
	AllocateFaces(P_FarfieldFaceNum,P_FarfieldFaceA,P_FarfieldVec);
	AllocateFaces(V_InletFaceNum,V_InletFaceA,V_InletVec);
	AllocateFaces(BoundFaceNum,BoundFaceA,BoundVec);
}
void ShadowCellCornerConstruct()
{
	double Halfdx = MinL/2.0;
	double Halfdy = MinL/2.0;
	PeriodicShadowC_NE = new Cell_2D();
	PeriodicShadowC_NW = new Cell_2D();
	PeriodicShadowC_SE = new Cell_2D();
	PeriodicShadowC_SW = new Cell_2D();
	for(int n = 0;n < PeriodicFaceNum;++n)
	{
		if
		(
			(fabs(PeriodicFaceA[n]->yf - Y_Beg - Halfdy) < infinitesimal)
			&&
			(X_Beg == PeriodicFaceA[n]->xf)
		)
		{
			PeriodicShadowC_NE->ShadowC = PeriodicFaceA[n]->lhsCell;
			PeriodicShadowC_NE->xc = PeriodicShadowC_NE->ShadowC->xc + Lx;
			PeriodicShadowC_NE->yc = PeriodicShadowC_NE->ShadowC->yc + Ly;
			PeriodicShadowC_NE->volume = PeriodicShadowC_NE->ShadowC->volume;
			PeriodicShadowC_NE->celltype = PeriodicShadowC_NE->ShadowC->celltype;
		}
		else if
		(
			(fabs(PeriodicFaceA[n]->yf - Y_Beg - Halfdy) < infinitesimal)
			&&
			(X_End == PeriodicFaceA[n]->xf)
		)
		{
			PeriodicShadowC_NW->ShadowC = PeriodicFaceA[n]->lhsCell;
			PeriodicShadowC_NW->xc = PeriodicShadowC_NW->ShadowC->xc - Lx;
			PeriodicShadowC_NW->yc = PeriodicShadowC_NW->ShadowC->yc + Ly;
			PeriodicShadowC_NW->volume = PeriodicShadowC_NW->ShadowC->volume;
			PeriodicShadowC_NW->celltype = PeriodicShadowC_NW->ShadowC->celltype;
		}
		else if
		(
			(fabs(Y_End - Halfdy - PeriodicFaceA[n]->yf) < infinitesimal)
			&&
			(X_Beg == PeriodicFaceA[n]->xf)
		)
		{
			PeriodicShadowC_SE->ShadowC = PeriodicFaceA[n]->lhsCell;
			PeriodicShadowC_SE->xc = PeriodicShadowC_SE->ShadowC->xc + Lx;
			PeriodicShadowC_SE->yc = PeriodicShadowC_SE->ShadowC->yc - Ly;
			PeriodicShadowC_SE->volume = PeriodicShadowC_SE->ShadowC->volume;
			PeriodicShadowC_SE->celltype = PeriodicShadowC_SE->ShadowC->celltype;
		}
		else if
		(
			(fabs(Y_End - Halfdy - PeriodicFaceA[n]->yf) < infinitesimal)
			&&
			(X_End == PeriodicFaceA[n]->xf)
		)
		{
			PeriodicShadowC_SW->ShadowC = PeriodicFaceA[n]->lhsCell;
			PeriodicShadowC_SW->xc = PeriodicShadowC_SW->ShadowC->xc - Lx;
			PeriodicShadowC_SW->yc = PeriodicShadowC_SW->ShadowC->yc - Ly;
			PeriodicShadowC_SW->volume = PeriodicShadowC_SW->ShadowC->volume;
			PeriodicShadowC_SW->celltype = PeriodicShadowC_SW->ShadowC->celltype;
		}
	}
}
void ShadowCPeriodicConstruct(int PeriodicFaceNum)
{
	if(0 == PeriodicFaceNum) return;
	PeriodicShadowCA = new Cell_2D[PeriodicFaceNum];
	int k = 0;
	for(int i = 0;i != Faces;++i)
	{
		if(12 == FaceArray[i].bc_type || 8 == FaceArray[i].bc_type)
		{
			FaceArray[i].rhsCell = PeriodicShadowCA + k;
			PeriodicShadowCA[k].Face_C[0] =  FaceArray + i;
			PeriodicShadowCA[k].Cell_C[0] =  FaceArray[i].lhsCell;
			PeriodicShadowCA[k].ShadowC = FaceArray[i].shadowF->lhsCell;
			PeriodicShadowCA[k].xc = PeriodicShadowCA[k].ShadowC->xc;
			PeriodicShadowCA[k].yc = PeriodicShadowCA[k].ShadowC->yc;
			PeriodicShadowCA[k].volume = PeriodicShadowCA[k].ShadowC->volume;
			PeriodicShadowCA[k].celltype = PeriodicShadowCA[k].ShadowC->celltype;
			if(X_Beg == FaceArray[i].xf)
			{
				PeriodicShadowCA[k].xc -= Lx;
			}
			else if(X_End == FaceArray[i].xf)
			{
				PeriodicShadowCA[k].xc += Lx;
			}
			else if(Y_Beg == FaceArray[i].yf)
			{
				PeriodicShadowCA[k].yc -= Ly;
			}
			else if(Y_End == FaceArray[i].yf)
			{
				PeriodicShadowCA[k].yc += Ly;
			}
			else
			{
				cout <<"Construct Shadow Cell Failed"<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
			}
			++k;
		}
	}
	#ifdef _CARTESIAN_MESH_FLIP
	ShadowCellCornerConstruct();
	#endif
}
void ShadowCBoundConstruct(int BoundFaceNum,Face_2D** &ptrFaceA,Cell_2D* &ptrShadowCA)
{
	if(0 == BoundFaceNum) return;
	ptrShadowCA = new Cell_2D[BoundFaceNum]();
	for(int i = 0;i != BoundFaceNum;++i)
	{
		ptrFaceA[i]->rhsCell = ptrShadowCA + i;
		ptrShadowCA[i].Face_C[0] = ptrFaceA[i];
		//PushCell(*ptrFaceA[i]->lhsCell,ptrShadowCA + i);
		ptrShadowCA[i].ShadowC = ptrFaceA[i]->lhsCell;
		ptrShadowCA[i].Cell_C[0] = ptrFaceA[i]->lhsCell;
		if(3 == ptrFaceA[i]->bc_type)
		{
			ptrShadowCA[i].xc = ptrFaceA[i]->xf;
			ptrShadowCA[i].yc = ptrFaceA[i]->yf;
		}
		else
		{
			ptrShadowCA[i].xc = 2.0*ptrFaceA[i]->xf - ptrShadowCA[i].Cell_C[0]->xc;
			ptrShadowCA[i].yc = 2.0*ptrFaceA[i]->yf - ptrShadowCA[i].Cell_C[0]->yc;	
		}
		ptrShadowCA[i].celltype = ptrFaceA[i]->lhsCell->celltype;
		ptrShadowCA[i].volume = ptrFaceA[i]->lhsCell->volume;
	}
}
int ShadowCellConstruct()
{
	ShadowCPeriodicConstruct(PeriodicFaceNum);
	ShadowCBoundConstruct(WallFaceNum,WallFaceA,WallShadowCA);
	ShadowCBoundConstruct(P_InletFaceNum,P_InletFaceA,P_InletShadowCA);
	ShadowCBoundConstruct(P_OutletFaceNum,P_OutletFaceA,P_OutletShadowCA);
	ShadowCBoundConstruct(SymmetryFaceNum,SymmetryFaceA,SymmetryShadowCA);
	ShadowCBoundConstruct(P_FarfieldFaceNum,P_FarfieldFaceA,P_FarfieldShadowCA);
	ShadowCBoundConstruct(V_InletFaceNum,V_InletFaceA,V_InletShadowCA);
	return 0;
}
void ShadowCellCheck(int BoundFaceNum,const Cell_2D* const &BoundShadowCA,const string& s)
{
	for(int i = 0;i < BoundFaceNum;++i)
	{
		if(nullptr == BoundShadowCA[i].ShadowC)
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" ShadowC != nullptr"<<endl;
			getchar();
		}
		if(nullptr == BoundShadowCA[i].Cell_C[0])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" Cell_C[0] == nullptr"<<endl;
			getchar();
		}
		if(nullptr == BoundShadowCA[i].Face_C[0])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" Face_C[0] == nullptr"<<endl;
			getchar();
		}
	}
}
void sortFacesInThisCell(Cell_2D &cell)
{
	Face_2D *Facetmp[4] = {nullptr,nullptr,nullptr,nullptr};
	for(int i = 0;i < cell.celltype;++i)
	{
		if(fabs(cell.Face_C[i]->xf - cell.xc) > infinitesimal)
		{
			if(cell.Face_C[i]->xf > cell.xc)
			{
				Facetmp[0] = cell.Face_C[i];
			}
			else
			{
				Facetmp[2] = cell.Face_C[i];
			}
		}
		else
		{
			if(cell.Face_C[i]->yf > cell.yc)
			{
				Facetmp[1] = cell.Face_C[i];
			}
			else
			{
				Facetmp[3] = cell.Face_C[i];
			}
		}
	}
	for(int i = 0;i < cell.celltype;++i)
	{
		cell.Face_C[i] = Facetmp[i];
	}
}
void NeighbourCellConstruct()
{
	for(int n = 0;n < Cells;++n)
	for(int i = 0;i < CellArray[n].celltype;++i)
	{
	#ifdef _CARTESIAN_MESH_FLIP
		sortFacesInThisCell(CellArray[n]);
	#endif
	CellArray[n].Cell_C[i] = ((CellArray[n].Face_C[i] -> lhsCell ==  (CellArray + n)) ? 
					CellArray[n].Face_C[i] -> rhsCell : CellArray[n].Face_C[i] -> lhsCell);
	CellArray[n].signFlux[i] = ((CellArray[n].Face_C[i] -> lhsCell ==  (CellArray + n)) ? -1 : 1);
	}
}
void DiagonalCellConstruct()
{
	int ne = 0,nw = 0,se = 0,sw = 0;
	for(int n = 0;n < Cells;++n)
	{
	//---------------------------------Diagonal 0 3------------------------
		if(nullptr == CellArray[n].Cell_C[0]->ShadowC)
		{
			CellArray[n].Cell_Diag[0] = CellArray[n].Cell_C[0]->Cell_C[1];
			CellArray[n].Cell_Diag[3] = CellArray[n].Cell_C[0]->Cell_C[3];
		}
		else
		{
			if(nullptr == CellArray[n].Cell_C[1]->ShadowC)
			{
				CellArray[n].Cell_Diag[0] = CellArray[n].Cell_C[1]->Cell_C[0];
			}
			else
			{
				CellArray[n].Cell_Diag[0] = PeriodicShadowC_NE;
				++ne;
			}
			if(nullptr == CellArray[n].Cell_C[3]->ShadowC)
			{
				CellArray[n].Cell_Diag[3] = CellArray[n].Cell_C[3]->Cell_C[0];
			}
			else
			{
				CellArray[n].Cell_Diag[3] = PeriodicShadowC_SE;
				++se;
			}
		}
	//---------------------------------Diagonal 1 2------------------------
		if(nullptr == CellArray[n].Cell_C[2]->ShadowC)
		{
			CellArray[n].Cell_Diag[1] = CellArray[n].Cell_C[2]->Cell_C[1];
			CellArray[n].Cell_Diag[2] = CellArray[n].Cell_C[2]->Cell_C[3];
		}
		else
		{
			if(nullptr == CellArray[n].Cell_C[1]->ShadowC)
			{
				CellArray[n].Cell_Diag[1] = CellArray[n].Cell_C[1]->Cell_C[2];
			}
			else
			{
				CellArray[n].Cell_Diag[1] = PeriodicShadowC_NW;
				++nw;
			}
			if(nullptr == CellArray[n].Cell_C[3]->ShadowC)
			{
				CellArray[n].Cell_Diag[2] = CellArray[n].Cell_C[3]->Cell_C[2];
			}
			else
			{
				CellArray[n].Cell_Diag[2] = PeriodicShadowC_SW;
				++sw;
			}
		}
	//---------------------------------Diagonal 2------------------------
	}
}
void LSCellMatrix(Cell_2D* const &Center,int k,
					const double &neighbour_xc, double const &neighbour_yc)
{
	double dx = neighbour_xc - Center->xc, dy = neighbour_yc - Center->yc;
//
	SetZero(dx);
	SetZero(dy);
//-----------------LeastSquare Left Matrix----------------------
	double DistanceP2P = dx*dx + dy*dy;
	Center->LS_M[0][0] += dx*dx/DistanceP2P;
	Center->LS_M[0][1] += dx*dy/DistanceP2P;
	Center->LS_M[1][0] += dy*dx/DistanceP2P;
	Center->LS_M[1][1] += dy*dy/DistanceP2P;
//-----------------LeastSquare Right Matrix---------------------
	Center->wdx_C[k] = dx/DistanceP2P;
	Center->wdy_C[k] = dy/DistanceP2P;
}

void InverseMatrix_2_2(double (&LS_M)[2][2])
{
	double a[2][2],A;
	A = LS_M[0][0] * LS_M[1][1] - LS_M[0][1] * LS_M[1][0];
	if(0.0 == A)
	{
		cout <<"Singular Maxtrix : " <<__FILE__<<"  "<<__LINE__<<"  "<<__func__<<endl;
		getchar();
		return;
	}
	a[0][0] = LS_M[1][1]/A;
	a[1][1] = LS_M[0][0]/A;
	a[0][1] = -LS_M[0][1]/A;
	a[1][0] = -LS_M[1][0]/A;
	for(int i = 0;i < 2;++i)
		for(int j = 0;j < 2;++j)
			LS_M[i][j] = a[i][j];
}

void Grad_LSMatrix()
{
	Cell_2D* neighbour = nullptr;
	Cell_2D* center = nullptr;
	for(int i = 0;i != Cells;++i)
	{
		center = CellArray + i;
		for(int k = 0;k != center->celltype;++k)
		{
			neighbour = center->Cell_C[k];
			if(neighbour != nullptr)
			{
				LSCellMatrix(center,k,neighbour->xc,neighbour->yc);
			}
			else
			{
				cout << "CellArray : " << i <<"neighbour cell invalid : "<<endl;
				getchar();
			}
		}
		InverseMatrix_2_2(center->LS_M);
	}
	cout <<"LeastSquare Matrix Construction Done" << endl;
}
void setXiDotdS()
{
	for(int n = 0;n < Faces;++n)
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		FaceArray[n].xi_n_dS[i][j] = FaceArray[n].Area*((xi_u[QuIndex])*FaceArray[n].Vx + xi_v[j]*FaceArray[n].Vy);
	}
}
int MeshCheck()
{
	cout <<"Face Checking..." <<endl;
	for(int i = 0;i != Faces;++i)
	{
		/*if(FaceArray[i].bc_type != 12 &&  FaceArray[i].bc_type != 8)
		cout << "FaceArray : " << i <<"   bc_type : " << FaceArray[i].bc_type <<std::hex
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[0] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[1] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].lhsCell , CellArray)
			 <<"  "<< MeshIndex(FaceArray[i].rhsCell , CellArray)<< std::dec <<endl;
		else
		cout << "FaceArray : " << i <<"   bc_type : " << FaceArray[i].bc_type <<std::hex
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[0] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[1] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].lhsCell , CellArray)
			 <<"  "<< MeshIndex(FaceArray[i].rhsCell->ShadowC , CellArray)<< std::dec <<endl;*/	
		#ifdef _CARTESIAN_MESH_FLIP
		if(FaceArray[i].Vx * FaceArray[i].Vy != 0.0)
		{
			cout << "FaceArray : " << i <<" Vx : " << FaceArray[i].Vx
				 <<" Vx : " << FaceArray[i].Vy <<" Area: "
				 <<std::setiosflags(ios::scientific) <<std::setprecision(16)
				 <<FaceArray[i].Area
				 <<std::resetiosflags(ios::scientific) <<std::setprecision(6)
				 << endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
		}
		if(FaceArray[i].Area != MinL)
		{
			cout << "FaceArray : " << i <<" Vx : " << FaceArray[i].Vx
				 <<" Vx : " << FaceArray[i].Vy <<" Area: "
				 <<std::setiosflags(ios::scientific) <<std::setprecision(16)
				 <<FaceArray[i].Area
				 <<std::resetiosflags(ios::scientific) <<std::setprecision(6)
				 << endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
		}
		#endif
		if(2 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"Cell = -1" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"construct shadow faces for interior faces"<< endl;
				getchar();
			}
		}
		else if(3 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"lhsCell == nullptr || rhsCell == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for wall boundaries"<< endl;
				getchar();
			}
		}
		else if(8 == FaceArray[i].bc_type || 12 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"Cell = -1" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"Failed to construct shadowF faces for 8 & 12"<< endl;
				getchar();
			}
		}
		else if(4 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"lhsCell == nullptr || rhsCell == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Pressure Inlet boundaries"<< endl;
				getchar();
			}
		}
		else if(5 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"lhsCell == nullptr || rhsCell == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Pressure Outlet boundaries"<< endl;
				getchar();
			}
		}
		else if(7 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"lhsCell == nullptr || rhsCell == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Symmetry boundaries"<< endl;
				getchar();
			}
		}
		else if(9 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"lhsCell == nullptr || rhsCell == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Pressure Farfield boundaries"<< endl;
				getchar();
			}
		}
		else if(10 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].lhsCell == nullptr || FaceArray[i].rhsCell == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"lhsCell == nullptr || rhsCell == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Velocity Inlet boundaries"<< endl;
				getchar();
			}
		}
		else
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "FaceArray : " << i <<"unknown bc_type : " << FaceArray[i].bc_type << endl;
			getchar();
		}
	}
	cout <<"Face Check Done" <<endl;
	cout <<"Cell Checking..." <<endl;
	for(int i = 0;i < Cells;++i)
	{
		for(int iF = 0;iF < CellArray[i].celltype;++iF)
		{
			if(nullptr == CellArray[i].Cell_C[iF])
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"Face Addr : "<<CellArray[i].Face_C[iF]<<'\n'
					<<"Face lhsCell : "<<CellArray[i].Face_C[iF]->lhsCell<<'\n'
					<<"Face rhsCell : "<<CellArray[i].Face_C[iF]->rhsCell<<'\n'
					<<"Cell Addr : "<<CellArray + i<<'\n';
			}
			if(CellArray[i].Face_C[iF]->lhsCell == CellArray + i)
			{
				if(CellArray[i].Face_C[iF]->rhsCell != CellArray[i].Cell_C[iF])
				{
					_PRINT_ERROR_MSG_FLIP
					cout <<"Construct neighbour cells failed : "<<"Face_C[iF] != Cell_C[iF]"<<endl
						 <<"CellArray : " << i <<endl;
					cout<<"Face_C[iF]->rhsCell : "<<CellArray[i].Face_C[iF]->rhsCell<<endl
						<<"Cell_C[iF]                      : "<<CellArray[i].Cell_C[iF]<<endl;
					getchar();
				}
			}
			else if(CellArray[i].Face_C[iF]->rhsCell == CellArray + i)
			{
				if(CellArray[i].Face_C[iF]->lhsCell != CellArray[i].Cell_C[iF])
				{
					_PRINT_ERROR_MSG_FLIP
					cout <<"Construct neighbour cells failed : "<<"Face_C[iF] != Cell_C[iF]"<<endl
						 <<"CellArray : " << i <<endl;
					cout<<"Face_C[iF]->lhsCell : "<<CellArray[i].Face_C[iF]->lhsCell<<endl
						<<"Cell_C[iF]                      : "<<CellArray[i].Cell_C[iF]<<endl;
					getchar();
				}
			}
			else
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Construct neighbour cells failed"<<endl<<"CellArray : " << i <<endl;
				getchar();
			}
		}
		if(nullptr != CellArray[i].ShadowC)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"CellArray : " << i <<endl;
			getchar();
		}
		if(3 == CellArray[i].celltype)
		{
			for(int k = 0;k != NumPerCell;++k)
			if(nullptr == CellArray[i].NodeX_C[k] || nullptr == CellArray[i].NodeY_C[k])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Node Index : " << k <<" = nullptr"<<endl;
				getchar();
			}

			if((CellArray[i].NodeX_C[2] != CellArray[i].NodeX_C[3]) ||
			   (CellArray[i].NodeY_C[2] != CellArray[i].NodeY_C[3]))
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Addrs of last two Nodes don't equal"<<endl;
				getchar();
			}

			for(int k = 0;k != NumPerCell - 1;++k)
			if(nullptr == CellArray[i].Cell_C[k])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Cell Index : " << k <<" = nullptr"<<endl;
				getchar();
			}

			if(nullptr != CellArray[i].Cell_C[NumPerCell - 1])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Cell Index : " << NumPerCell - 1 <<" != nullptr"<<endl;
				getchar();
			}
		}
		else if(4 == CellArray[i].celltype)
		{
			for(int k = 0;k != NumPerCell;++k)
			if(nullptr == CellArray[i].NodeX_C[k] || nullptr == CellArray[i].NodeY_C[k])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Quadrilateral Cell "<<"CellArray : " << i 
					 <<"  Node Index : " << k <<" == nullptr"<<endl;
				getchar();
			}
			for(int k = 0;k != NumPerCell;++k)
			if(nullptr == CellArray[i].Cell_C[k])
			{
				cout <<"Quadrilateral Cell "<<"CellArray : " << i 
					 <<"  Cell Index : " << k <<" == nullptr"<<endl;
				getchar();
			}
		}
		else
		{
			cout << "CellArray : " << i <<" unknown celltype : "
			     << CellArray[i].celltype << endl;
			getchar();
		}
		#ifdef _CARTESIAN_MESH_FLIP
		if(MinL*MinL != CellArray[i].volume)
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "CellArray : " << i <<" unmatch cell volume : "
				 <<std::setiosflags(ios::scientific) <<std::setprecision(16)
			     << MinL*MinL <<"  "<< CellArray[i].volume <<"    "<<MinL*MinL - CellArray[i].volume
			     <<endl;
			getchar();
		}
		#endif
	}
	for(int n = 0;n < PeriodicFaceNum;++n)
	{
		const string &s = "PeriodicShadowCA";
		if(nullptr == PeriodicShadowCA[n].ShadowC)
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<n<<" ShadowC != nullptr"<<endl;
			getchar();
		}
		for(int k = 1;k < NumPerCell;++k)		
		if(nullptr != PeriodicShadowCA[n].Cell_C[k])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<n<<" Cell_C[k] == nullptr"<<endl;
			getchar();
		}
		for(int k = 1;k < NumPerCell;++k)
		if(nullptr != PeriodicShadowCA[n].Face_C[k])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<n<<" Face_C[k] == nullptr"<<endl;
			getchar();
		}
	}
	ShadowCellCheck(WallFaceNum,WallShadowCA,"WallShadowCA");
	ShadowCellCheck(P_InletFaceNum,P_InletShadowCA,"P_InletShadowCA");
	ShadowCellCheck(P_OutletFaceNum,P_OutletShadowCA,"P_OutletShadowCA");
	ShadowCellCheck(SymmetryFaceNum,SymmetryShadowCA,"SymmetryShadowCA");
	ShadowCellCheck(P_FarfieldFaceNum,P_FarfieldShadowCA,"P_FarfieldShadowCA");
	ShadowCellCheck(V_InletFaceNum,V_InletShadowCA,"V_InletShadowCA");
	cout <<"Cell Check Done" <<endl;
	return 0;
}
