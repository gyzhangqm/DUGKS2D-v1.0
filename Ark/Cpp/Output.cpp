#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <string>
#include "DUGKSDeclaration.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::map;
using std::ios;
using std::setprecision;

int const Out_precision = 15;

int const IC = 198, JC = 197;

extern void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p);

extern void TaylorCouetteAnalyticalSolution(double x,double y,double &u_A);

extern void AnalyticalForceDrivenTG(double x,double y,double &u_A, double &v_A,double &p_A);

void OutputCase()
{
	ofstream OutFile_Case("../FlowField/Case.ark");
	if(!OutFile_Case)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"Case.ark open failed."<<endl;
		getchar();
		return;
	}
	OutFile_Case <<_MESHTYPE_ARK<<"    =    Mesh Type"<<endl
				 <<_FLUX_SCHEME_ARK<<"    =    Flux Scheme"<<endl
				 <<_BC_ARK<<"    =    Boundary Condition"<<endl
				 <<_QMODEL_ARK<<"    =    velocity model"<<endl
				 <<_MESHFILE_NAME_ARK<<"    =    Mesh File"<<endl
				 <<"#----------------------------------------"<<endl
				 <<Omega0<<"    =    Omega0//sutherland power"<<endl
				 <<Pr<<"    =    Prandtl"<<endl
				 <<nK<<"    =    nK//internal degrees of freedom"<<endl
				 <<Cv<<"    =    specific heat of constant volume"<<endl
				 <<Gamma<<"    =    Gamma//specific heat ratio"<<endl
				 <<"#----------------------------------------"<<endl
				 <<DV_Qu<<"    =    Qu"<<endl
				 <<DV_Qv<<"    =    Qv"<<endl
				 <<NL<<"    =    NL"<<endl
				 <<ChLength<<"    =    Character Length"<<endl
				 <<MinL<<"    =    MinL"<<endl
				 <<Kn<<"    =    Knudsen"<<endl
				 <<Ma<<"    =    Mach"<<endl
				 <<Re<<"    =    Re"<<endl
				 <<U0<<"    =    U0"<<endl
				 <<R0<<"    =    R0"<<endl
				 <<T0<<"    =    T0"<<endl
				 <<Mu0<<"    =    dynamic viscosity"<<endl
				 <<Lambda0<<"    =    Lambda0"<<endl
				 <<CFL<<"    =    CFL"<<endl
				 <<RESIDUAL<<"    =    RESIDUAL"<<endl
				 <<dt<<"    =    dt"<<endl
				 <<Mu0/(Rho0*R0*T0)<<"    =    tau0"<<endl
				 <<dt*Rho0*R0*T0/Mu0<<"    =    dt/tau0"<<endl
				 <<MaSpan<<"    =    MaSpan"<<endl
				 <<Eta<<"    =    Eta// T of D2V16";
//
	cout 		 <<_MESHTYPE_ARK<<"    =    Mesh Type"<<endl
				 <<_FLUX_SCHEME_ARK<<"    =    Flux Scheme"<<endl
				 <<_BC_ARK<<"    =    Boundary Condition"<<endl
				 <<_QMODEL_ARK<<"    =    velocity model"<<endl
				 <<_MESHFILE_NAME_ARK<<"    =    Mesh File"<<endl
				 <<"#----------------------------------------"<<endl
				 <<NL<<"    =    NL"<<endl
				 <<ChLength<<"    =    Character Length"<<endl
				 <<MinL<<"    =    MinL"<<endl
				 <<Kn<<"    =    Knudsen"<<endl
				 <<Ma<<"    =    Mach"<<endl
				 <<Re<<"    =    Re"<<endl
				 <<U0<<"    =    U0"<<endl
				 <<R0<<"    =    R0"<<endl
				 <<T0<<"    =    T0"<<endl
				 <<Mu0<<"    =    dynamic viscosity"<<endl
				 <<Lambda0<<"    =    Lambda0"<<endl
				 <<CFL<<"    =    CFL"<<endl
				 <<RESIDUAL<<"    =    RESIDUAL"<<endl
				 <<dt<<"    =    dt"<<endl
				 <<Mu0/(Rho0*R0*T0)<<"    =    tau0"<<endl
				 <<dt*Rho0*R0*T0/Mu0<<"    =    dt/tau0"<<endl
				 <<MaSpan<<"    =    MaSpan"<<endl
				 <<Eta<<"    =    Eta// T of D2V16";
	OutFile_Case.close();
}
void Output_Convergence()
{
	ostringstream oss_L2;
	oss_L2 <<"../Convergence/L2_uvp_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat"; 
	ofstream OutFile_L2(oss_L2.str().c_str());
	if(!OutFile_L2)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_L2 file open failed"<<endl;
		getchar();
	}
//---------------------------------------------------------------------
	ostringstream oss_SumRho;
	oss_SumRho << "../Convergence/SumRho_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_SumRho(oss_SumRho.str().c_str());
	if(!OutFile_SumRho)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_SumRho file open failed"<<endl;
		getchar();
	}
//---------------------------------------------------------------------
	ostringstream oss_Residual;
	oss_Residual <<"../Convergence/Residual_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_Residual(oss_Residual.str().c_str());
	if(!OutFile_Residual)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_Residual open failed"<<endl;
		getchar();
	}
}
void TaylorGreen_L2Norm(double const &t,double &L2_uv, double &L2_p)
{
	double u_A, v_A, p_A, du, dv, dp; 
	double  Sumdudv = 0.0,Sumdp = 0.0,Sumuv_A = 0.0,Sump_A = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A,v_A,p_A);
		du = CellArray[i].U - u_A;
		dv = CellArray[i].V - v_A;
		dp = CellArray[i].p - p_A;
		Sumdudv += du*du + dv*dv;
		Sumuv_A += u_A*u_A + v_A*v_A; 
		Sumdp   += dp*dp;
		Sump_A  += p_A*p_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
	L2_p  = sqrt(Sumdp/Sump_A);
}
void TaylorCouette_L2Norm(double const &t,double &L2_uv)
{
	double u_A, du;
	double  Sumdudv = 0.0,Sumuv_A = 0.0;
	for(int i =0;i < Cells;++i)
	{
		TaylorCouetteAnalyticalSolution(CellArray[i].xc,CellArray[i].yc,u_A);
		du = sqrt(CellArray[i].U*CellArray[i].U + CellArray[i].V*CellArray[i].V) - u_A;
		Sumdudv += du*du;
		Sumuv_A += u_A*u_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
}
void ForceDrivenTaylorGreen_L2Norm(double &L2_uv, double &L2_p)
{
	double u_A, v_A, p_A, du, dv, dp;
	double  Sumdudv = 0.0,Sumdp = 0.0,Sumuv_A = 0.0,Sump_A = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		AnalyticalForceDrivenTG(CellArray[i].xc,CellArray[i].yc,u_A,v_A,p_A);
		du = CellArray[i].U - u_A;
		dv = CellArray[i].V - v_A;
		dp = CellArray[i].p - p_A;
		Sumdudv += du*du + dv*dv;
		Sumuv_A += u_A*u_A + v_A*v_A; 
		Sumdp   += dp*dp;
		Sump_A  += p_A*p_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
	L2_p  = sqrt(Sumdp/Sump_A);
}
void Output_L2Norm(double const &t,double &L2_uv, double &L2_p)
{	
	//ForceDrivenTaylorGreen_L2Norm(L2_uv,L2_p);
	TaylorCouette_L2Norm(t,L2_uv);
	ostringstream oss_L2;
	oss_L2 <<"../Convergence/L2_uvp_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_L2(oss_L2.str().c_str(),ofstream::app);
	if(!OutFile_L2)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_L2 file open failed"<<endl;
		getchar();
	}
	OutFile_L2 << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_L2 << t <<"    "<<L2_uv<<"    "<<L2_p<<'\n';
	OutFile_L2.close();
}
void Output_SumRho(double t)
{
	ostringstream oss_SumRho;
	oss_SumRho << "../Convergence/SumRho_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_SumRho(oss_SumRho.str().c_str(),ofstream::app);
	if(!OutFile_SumRho)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_SumRho file open failed"<<endl;
		getchar();
	}
	OutFile_SumRho << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_SumRho <<  t <<"    "<<SumRho<<endl;
	OutFile_SumRho.close();
}
void Output_Residual(double t,double Residual)
{
	ostringstream oss_Residual;
	oss_Residual <<"../Convergence/Residual_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_Residual(oss_Residual.str().c_str(),ofstream::app);
	if(!OutFile_Residual)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_Residual open failed"<<endl;
		getchar();
	}
	OutFile_Residual << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_Residual << t <<"    "<<Residual<<endl;
	OutFile_Residual.close();
}
/*void writeZone(const string &xyz ,const string &uvp,
				const map<double,double> &profile,double const t,ofstream &OutFile_profile)
{
	OutFile_profile <<"VARIABLES ="<< xyz <<","<< uvp <<endl;
	OutFile_profile <<"ZONE T = Time" << t <<"_"<<uvp<<endl;
	OutFile_profile <<"I = "<<profile.size()<<",F = POINT"<<endl;
	for(auto p_It = profile.cbegin();p_It != profile.cend();++p_It)
		OutFile_profile << p_It->first <<"    "<<p_It->second<<endl;
}
void writeProfile(const string &xyz ,const string &uvp,double const t,
					const map<double,double> &profile,const map<double,double> &profile_A)
{
	ostringstream oss_profile;
	string uvp_A = uvp + "_A";
	oss_profile <<"../FlowField/"<< uvp <<"_MeshCar" <<NL<<"-"<<NL
				<<"_nu"<<nu<<"_RT"<<RT<<"_dt"<<CFL<<".dat";
	ofstream OutFile_profile(oss_profile.str().c_str());
	if(!OutFile_profile)
	{
		cout << uvp <<" profile open failed"<<endl;
		_PRINT_ERROR_MSG_FLIP
		getchar();
		return;
	}
	if(profile.size() != profile_A.size())
	{
		cout << uvp <<" : map size doesn't equal" <<endl;
		getchar();
		return;
	}
	OutFile_profile <<setiosflags(ios::scientific)<<setprecision(12);
	writeZone(xyz,uvp,profile,t,OutFile_profile);
	writeZone(xyz,uvp_A,profile_A,t,OutFile_profile);
	OutFile_profile.close();	
}*/
void Output_Flowfield(double const &t,int step)
{
	ostringstream oss_FlowField;
	oss_FlowField <<"../FlowField/global/" << "Step" << step <<"Ma"<< Ma<<"_"
					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"_T"<<t<<".dat";
	ofstream OutFile_FlowField(oss_FlowField.str().c_str());
	if(!OutFile_FlowField)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_FlowField open failed" << endl; 
		getchar();
		return;
	}
	ostringstream VarName,VarLocation,ZoneName,dataNE;
	VarName << "VARIABLES = X,Y,<Greek>r</Greek>,U,V,p,T,q<sub>x</sub>,q<sub>y</sub>,\
	<Greek>t</Greek><sub>xx</sub>,<Greek>t</Greek><sub>xy</sub>,<Greek>t</Greek><sub>yy</sub>\n";
	VarLocation <<"VarLocation=([1-2]=NODAL,[3-12]=CellCentered)\n";
	ZoneName<<"ZONE T = Time" << t <<"_"<<"Mu"<<Mu0<<"\n";
	dataNE<<"Nodes="<<Nodes<<", Elements="<<Cells<<", ZONETYPE=FEQuadrilateral\n";
	string tecformat[5]={VarName.str().c_str(),
						ZoneName.str().c_str(),
						dataNE.str().c_str(),
						"DATAPACKING=BLOCK\n",
						VarLocation.str().c_str()};
	OutFile_FlowField << tecformat[0]<<tecformat[1]<<tecformat[2]<<tecformat[3]<<tecformat[4];
	OutFile_FlowField << setiosflags(ios::scientific) << setprecision(12);
//	
	/*double *u_A = new double[Cells];
	double *v_A = new double[Cells];
	double *p_A = new double[Cells];
	for(int i = 0;i != Cells;++i)
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A[i], v_A[i], p_A[i]);
	*/
//--------------------------------------------NodeX-------------------------------------------
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_FlowField << NodeX[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	OutFile_FlowField << endl;
//--------------------------------------------NodeY-------------------------------------------
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_FlowField << NodeY[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	OutFile_FlowField << endl;
//--------------------------------------------rho-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].Rho <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------u-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].U <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------v-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].V <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------p-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].p <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].T <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].qx <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].qy <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].shearTau[0][0] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].shearTau[0][1] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].shearTau[1][1] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
/*//--------------------------------------------u_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << u_A[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------v_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << v_A[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------p_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << p_A[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}*/
//--------------------------------------------relation----------------------------------------	
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << MeshIndex(CellArray[i].NodeX_C[0] , NodeX)<< " "
					 << MeshIndex(CellArray[i].NodeX_C[1] , NodeX)<< " "
					 << MeshIndex(CellArray[i].NodeX_C[2] , NodeX)<< " " 
					 << MeshIndex(CellArray[i].NodeX_C[3] , NodeX)<< endl;
	}
	OutFile_FlowField.close();
	/*delete []u_A;
	delete []v_A;
	delete []p_A;*/
}
//----------------------------------------------DEBUG----------------------------------------
void oss_XXX(ostringstream& oss_r,const string &folder,const string &suffix,double const &t)
{
	oss_r << "../FlowField/"<<folder<<"/Time" <<t<<"."<<suffix;
}
void FileOpen(ofstream &OutFile_XXX,ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str());
	if(!OutFile_XXX)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<< ("OutFile_" + s + " open failed") << endl; 
		getchar();
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void FileOpenAppend(ofstream &OutFile_XXX,ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str(),ios::app);
	if(!OutFile_XXX)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<< ("OutFile_" + s + " open failed") << endl; 
		getchar();
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void Output_UVP(double const &t)
{
	ostringstream oss_rho,oss_u,oss_v,oss_p,oss_uA,oss_vA,oss_pA;
	ofstream OutFile_rho,OutFile_u, OutFile_v, OutFile_p, OutFile_uA, OutFile_vA, OutFile_pA;
//--------------------------------------------------------------------------------
	double *u_A = new double[Cells];
	double *v_A = new double[Cells];
	double *p_A = new double[Cells];
	for(int i = 0;i != Cells;++i)
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A[i], v_A[i], p_A[i]);
//---------------------------------------------------------------------------------
	oss_XXX(oss_rho,"UVP","rho",t);
	oss_XXX(oss_u,"UVP","u",t);
	oss_XXX(oss_v,"UVP","v",t);
	oss_XXX(oss_p,"UVP","p",t);
	oss_XXX(oss_uA,"UVP","uA",t);
	oss_XXX(oss_vA,"UVP","vA",t);
	oss_XXX(oss_pA,"UVP","pA",t);
	FileOpen(OutFile_rho,oss_rho,"rho");
	FileOpen(OutFile_u,oss_u,"u");
	FileOpen(OutFile_v,oss_v,"v");
	FileOpen(OutFile_p,oss_p,"p");
	FileOpen(OutFile_uA,oss_uA,"uA");
	FileOpen(OutFile_vA,oss_vA,"vA");
	FileOpen(OutFile_pA,oss_pA,"pA");
//--------------------------------------------rho-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_rho << CellArray[i].Rho <<"\n";
	}
//--------------------------------------------u-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_u << CellArray[i].U <<"\n";
	}
//--------------------------------------------v-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_v << CellArray[i].V <<"\n";
	}
//--------------------------------------------p-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_p << CellArray[i].p <<"\n";
	}
//--------------------------------------------u_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_uA << u_A[i] <<"\n";
	}
//--------------------------------------------v_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_vA << v_A[i] <<"\n";
	}
//--------------------------------------------p_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_pA << p_A[i] <<"\n";
	}
	OutFile_rho.close();
	OutFile_u.close();
	OutFile_v.close();
	OutFile_p.close();
	OutFile_uA.close();
	OutFile_vA.close();
	OutFile_pA.close();
	delete []u_A;
	delete []v_A;
	delete []p_A;
}
void Output_fBP(double const &t,int ii,int jj)
{
	ostringstream oss_fBP;
	oss_fBP <<"../FlowField/fBP/" << "Time" << t <<"_Mu"<< Mu0
					<<"_MeshCar"<<NL<<"-"<<NL<<"_CFL"<<CFL<<"_fBP"<<ii<<"_"<<jj<<".dat";
	ofstream OutFile_fBP(oss_fBP.str().c_str());
	if(!OutFile_fBP)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_fBP open failed" << endl;
		getchar();
		return;
	}
	OutFile_fBP << setiosflags(ios::scientific) << setprecision(12);
	for(int n = 0;n != Cells;++n)
	{
		OutFile_fBP << CellArray[n].fBP[ii][jj] <<"\n";
	}
	OutFile_fBP.close();
}
void Output_fBh(Face_2D& face,double t)
{
	ostringstream oss_fBh;
	ofstream OutFile_fBh;
	oss_XXX(oss_fBh,"fBh","fBh",t);
	FileOpen(OutFile_fBh,oss_fBh,"fBh");
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_fBh<<face.fBh[i][j]<<'\n';
		// if(face.xi_n_dS[i][j] >= 0)
		// {
		// 	OutFile_fBh <<face.lhsCell->fBP[i][j]<<"    "
		// 				<<face.lhsCell->fEq[i][j]<<"    "
		// 				<<face.fBh[i][j] - face.lhsCell->fBP[i][j]<<"    "
		// 				<<face.lhsCell->fEq[i][j] - face.lhsCell->fBP[i][j]<<endl;
		// }
		// else
		// {
		// 	OutFile_fBh <<face.rhsCell->fBP[i][j]<<"    "
		// 				<<face.rhsCell->fEq[i][j]<<"    "
		// 				<<face.fBh[i][j] - face.rhsCell->fBP[i][j]<<"    "
		// 				<<face.rhsCell->fEq[i][j] - face.rhsCell->fBP[i][j]<<endl;
		// }
	}
	OutFile_fBh.close();
}
void Output_gBh(Face_2D& face,double t)
{
	ostringstream oss_gBh;
	ofstream OutFile_gBh;
	oss_XXX(oss_gBh,"gBh","gBh",t);
	FileOpen(OutFile_gBh,oss_gBh,"gBh");
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_gBh<<face.gBh[i][j]<<'\n';
		// if(face.xi_n_dS[i][j] >= 0)
		// {
		// 	OutFile_gBh <<face.lhsCell->gBP[i][j]<<"    "
		// 				<<face.lhsCell->gEq[i][j]<<"    "
		// 				<<face.gBh[i][j] - face.lhsCell->gBP[i][j]<<"    "
		// 				<<face.lhsCell->gEq[i][j] - face.lhsCell->gBP[i][j]<<endl;
		// }
		// else
		// {
		// 	OutFile_gBh <<face.rhsCell->gBP[i][j]<<"    "
		// 				<<face.rhsCell->gEq[i][j]<<"    "
		// 				<<face.gBh[i][j] - face.rhsCell->gBP[i][j]<<"    "
		// 				<<face.rhsCell->gEq[i][j] - face.rhsCell->gBP[i][j]<<endl;
		// }
	}
	OutFile_gBh.close();
}
void Output_fh(Face_2D& face,double t)
{
	ostringstream oss_fh;
	ofstream OutFile_fh;
	oss_XXX(oss_fh,"fh","fh",t);
	FileOpen(OutFile_fh,oss_fh,"fh");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_fh << face.fh[i][j] <<"    "
					   << face.fBh[i][j]<<"    "
					   << face.fEqh[i][j]<<'\n';
		}
	OutFile_fh <<"Rho_h : "<<face.Rho_h<<'\n'
			   <<"U_h : "<<face.U_h<<'\n'
			   <<"V_h : "<<face.V_h<<'\n'
			   <<"T_h: "<<face.T_h<<'\n'
			   <<"qx_h : "<<face.qx_h<<'\n'
			   <<"qy_h : "<<face.qy_h<<'\n'
			   <<"Vx : "<<face.Vx<<'\n'
			   <<"Vy : "<<face.Vy<<'\n'
			   <<"xc : "<<face.xf<<'\n'
			   <<"yc : "<<face.yf<<'\n';
	OutFile_fh.close();
}
void Output_fh_Append(Face_2D& face,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_fh;
	ofstream OutFile_fh;
	oss_XXX(oss_fh,"fh","fhAPP",dt);
	FileOpenAppend(OutFile_fh,oss_fh,"fh");
	OutFile_fh 	<< face.fh[I][J] <<"    "<<face.ah<<"    "<<face.bh<<'\n';
				//<< face.fBh[I][J]<<"    "
				//<< face.fEqh[I][J]<<'\n';
	OutFile_fh.close();				
}
void Output_gh_Append(Face_2D& face,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_gh;
	ofstream OutFile_gh;
	oss_XXX(oss_gh,"gh","ghAPP",dt);
	FileOpenAppend(OutFile_gh,oss_gh,"gh");
	OutFile_gh 	<< face.gh[I][J] <<"    "<<'\n';
				//<< face.fBh[I][J]<<"    "
				//<< face.fEqh[I][J]<<'\n';
	OutFile_gh.close();				
}
void Output_gh(Face_2D& face,double t)
{
	ostringstream oss_gh;
	ofstream OutFile_gh;
	oss_XXX(oss_gh,"gh","gh",t);
	FileOpen(OutFile_gh,oss_gh,"gh");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_gh << face.gh[i][j] <<"    "
					   << face.gBh[i][j]<<"    "
					   << face.gEqh[i][j]<<'\n';
		}
	OutFile_gh.close();
}
void Output_fT(Cell_2D &cell,double t)
{
	ostringstream oss_fT;
	ofstream OutFile_fT;
	oss_XXX(oss_fT,"fT","fT",t);
	FileOpen(OutFile_fT,oss_fT,"fT");
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_fT <<cell.fT[i][j]<<"    "
				   <<cell.fEq[i][j]<<"    "
//				   <<cell.fFlux[i][j]<<"    "
				   <<cell.fT[i][j]-cell.fEq[i][j]
				   <<'\n';
	}
	OutFile_fT <<"Rho : "<<cell.Rho<<'\n'
				<<"U : "<<cell.U<<'\n'
				<<"V : "<<cell.V<<'\n'
				<<"T: "<<cell.T<<'\n'
				<<"Lambda: "<<cell.Lambda<<'\n'
				<<"qx : "<<cell.qx<<'\n'
				<<"qy : "<<cell.qy<<'\n'
				<<"tauxx : "<<cell.shearTau[0][0]<<'\n'
				<<"tauxy : "<<cell.shearTau[0][1]<<'\n'
				<<"tauyy : "<<cell.shearTau[1][1]<<'\n'
				<<"xc : "<<cell.xc<<'\n'
				<<"yc : "<<cell.yc<<'\n';
	OutFile_fT.close(); 
}
void Output_gT(Cell_2D &cell,double t)
{
	ostringstream oss_gT;
	ofstream OutFile_gT;
	oss_XXX(oss_gT,"gT","gT",t);
	FileOpen(OutFile_gT,oss_gT,"gT");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_gT <<cell.gT[i][j]<<"    "
//					   <<cell.fEq[i][j]<<"    "
//					   <<cell.fFlux[i][j]<<"    "
//					   <<cell.gT[i][j]-cell.fEq[i][j]-cell.DtSlashVolume*cell.fFlux[i][j]
					   <<endl;
		}
	OutFile_gT <<"Rho : "<<cell.Rho<<'\n'
				<<"U : "<<cell.U<<'\n'
				<<"V : "<<cell.V<<'\n'
				<<"T: "<<cell.T<<'\n'
				<<"Lambda: "<<cell.Lambda<<'\n'
				<<"qx : "<<cell.qx<<'\n'
				<<"qy : "<<cell.qy<<'\n'
				<<"xc : "<<cell.xc<<'\n'
				<<"yc : "<<cell.yc<<'\n';
	OutFile_gT.close(); 
}
void Output_fT_Append(Cell_2D &cell,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_fT;
	ofstream OutFile_fT;
	oss_XXX(oss_fT,"fT","fTAPP",dt);
	FileOpenAppend(OutFile_fT,oss_fT,"fT");
	OutFile_fT     <<cell.fT[I][J]<<"    "//<<cell.aBP<<"    "<<cell.bBP
				<<cell.Cell_C[0]->fT[I][J]<<"    "
				<<cell.Cell_C[1]->fT[I][J]<<"    "
				<<cell.Cell_C[2]->fT[I][J]<<"    "
				<<cell.Cell_C[3]->fT[I][J]<<"    "
				<<'\n';
	OutFile_fT.close();
}
void Output_gT_Append(Cell_2D &cell,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_gT;
	ofstream OutFile_gT;
	oss_XXX(oss_gT,"gT","gTAPP",dt);
	FileOpenAppend(OutFile_gT,oss_gT,"gT");
	OutFile_gT  <<cell.gT[I][J]<<"    "
				// <<cell.gBP[I][J]<<"    "
				// <<cell.gEq[I][J]<<"    "
				<<cell.Cell_C[0]->gT[I][J]<<"    "
				<<cell.Cell_C[1]->gT[I][J]<<"    "
				<<cell.Cell_C[2]->gT[I][J]<<"    "
				<<cell.Cell_C[3]->gT[I][J]
				<<'\n';
	OutFile_gT.close();
}
// void Output_fFlux_Append(Cell_2D &cell,double dt)
// {
// 	int I = IC,J = JC;
// 	ostringstream oss_fFlux;
// 	ofstream OutFile_fFlux;
// 	oss_XXX(oss_fFlux,"fFlux","fFluxAPP",dt);
// 	FileOpenAppend(OutFile_fFlux,oss_fFlux,"fFlux");
// 	OutFile_fFlux  <<cell.fFlux[I][J]<<"    "
// 				    <<cell.Face_C[0]->fh[I][J]<<"    "
// 				    <<cell.Face_C[1]->fh[I][J]<<"    "
// 				    <<cell.Face_C[2]->fh[I][J]<<"    "
// 				    <<cell.Face_C[3]->fh[I][J]<<"    "
// 				    <<'\n';
// 	OutFile_fFlux.close();
// }
// void Output_gFlux_Append(Cell_2D &cell,double dt)
// {
// 	int I = IC,J = JC;
// 	ostringstream oss_gFlux;
// 	ofstream OutFile_gFlux;
// 	oss_XXX(oss_gFlux,"gFlux","gFluxAPP",dt);
// 	FileOpenAppend(OutFile_gFlux,oss_gFlux,"gFlux");
// 	OutFile_gFlux  <<cell.gFlux[I][J]<<"    "
// 				    <<cell.Face_C[0]->gh[I][J]<<"    "
// 				    <<cell.Face_C[1]->gh[I][J]<<"    "
// 				    <<cell.Face_C[2]->gh[I][J]<<"    "
// 				    <<cell.Face_C[3]->gh[I][J]<<"    "
// 				    <<'\n';
// 	OutFile_gFlux.close();
// }
void Output_phi_Bh(Face_2D &face,double t)
{
	ostringstream oss_fBh_L;
	ostringstream oss_fBh_R;
	ostringstream oss_gBh_L;
	ostringstream oss_gBh_R;
	ofstream OutFile_fBh_L;
	ofstream OutFile_fBh_R;
	ofstream OutFile_gBh_L;
	ofstream OutFile_gBh_R;
	oss_XXX(oss_fBh_L,"fBh","fBhL",t);
	oss_XXX(oss_fBh_R,"fBh","fBhR",t);
	oss_XXX(oss_gBh_L,"fBh","gBhL",t);
	oss_XXX(oss_gBh_R,"fBh","gBhR",t);
	FileOpen(OutFile_fBh_L,oss_fBh_L,"fBh_L");
	FileOpen(OutFile_fBh_R,oss_fBh_R,"fBh_R");
	FileOpen(OutFile_gBh_L,oss_gBh_L,"gBh_L");
	FileOpen(OutFile_gBh_R,oss_gBh_R,"gBh_R");
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		if(xi_u[QuIndex]>=0)
		{
			OutFile_fBh_R<<face.fBh[i][j]<<'\n';
			OutFile_gBh_R<<face.gBh[i][j]<<'\n';
		}
		else
		{
			OutFile_fBh_L<<face.fBh[i][j]<<'\n';
			OutFile_gBh_L<<face.gBh[i][j]<<'\n';
		}
	}
	OutFile_fBh_L.close();
	OutFile_fBh_R.close();
	OutFile_gBh_L.close();
	OutFile_gBh_R.close();
}
void Output_xcyc()
{
	ostringstream oss_xcyc;
	oss_xcyc <<"../FlowField/UVP/" << "Mu"<< Mu0
					<<"_MeshCar"<<NL<<"-"<<NL<<"_CFL"<<CFL<<"_xcyc.dat";
	ofstream OutFile_xcyc(oss_xcyc.str().c_str());
	if(!OutFile_xcyc)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_xcyc open failed" << endl; 
		getchar();
		return;
	}
	OutFile_xcyc << setiosflags(ios::scientific) << setprecision(16);
	for(int i = 0;i != Cells;++i)
	{
		OutFile_xcyc << CellArray[i].xc<<"  "<<CellArray[i].yc<<"\n";
	}
	OutFile_xcyc <<"--------------P_Inlet-----------------"<<'\n';
	for(int n = 0;n != P_InletFaceNum;++n)
		OutFile_xcyc << P_InletShadowCA[n].xc<<"  "<<P_InletShadowCA[n].yc<<"\n";
	OutFile_xcyc <<"--------------P_Outlet----------------"<<'\n';
	for(int n = 0;n != P_OutletFaceNum;++n)
		OutFile_xcyc << P_OutletShadowCA[n].xc<<"  "<<P_OutletShadowCA[n].yc<<"\n";
	for(int n = 0;n != PeriodicFaceNum;++n)
		OutFile_xcyc << PeriodicShadowCA[n].xc<<"  "<<PeriodicShadowCA[n].yc<<"\n";
	OutFile_xcyc.close();
}