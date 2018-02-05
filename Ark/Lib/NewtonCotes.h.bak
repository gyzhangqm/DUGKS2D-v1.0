#ifndef _ARK_NEWTONCOTES_H
#define _ARK_NEWTONCOTES_H
//
double CompositeTrapezoid(double** const &fEq,int const &Range,int const &j,double const &dx)
{
	double Sumu = 0.5*(fEq[0][j] + fEq[Range - 1][j]);
	for(int i = 1;i < Range - 1;++i)
	{
		Sumu += fEq[i][j];
	}
	return Sumu*dx;
}
double CompositeTrapezoid(double* const &fEq,int const &Range,double const &dx)
{
	double Sum = 0.5*(fEq[0] + fEq[Range - 1]);
	for(int i = 1;i < Range - 1;++i)
	{
		Sum += fEq[i];
	}
	return Sum*dx;
}
double CompositeTrapezoid(double* const &fEq,int const &Range,double* const &xi,double const &dx)
{
	double Sum_ru = 0.5*(fEq[0]*xi[0] + fEq[Range - 1]*xi[Range - 1]);
	for(int i = 1;i < Range - 1;++i)
	{
		Sum_ru += fEq[i]*xi[i];
	}
	return Sum_ru*dx;
}
double CompositeTrapezoid(double* const &fEq,int const &Range,
	double* const &xi_x,double* const &xi_y,double const &dx)
{
	double Sum = 0.5*(fEq[0]*xi_x[0]*xi_y[0] + fEq[Range - 1]*xi_x[Range - 1]*xi_y[Range - 1]);
	for(int i = 1;i < Range - 1;++i)
	{
		Sum += fEq[i]*xi_x[i]*xi_y[i];
	}
	return Sum*dx;
}
//-------------------------------------------------------------------------
double CompositeSimpson(double* const &fEq,int const &End,double const &dx)
{
	double sum = fEq[0] - fEq[End-1];
	double sumA = 0.0,sumB = 0.0;
	for(int i = 1;i < End-1;i+=2)
	{
		sumA += fEq[i];
		sumB += fEq[i+1];
	}
	sum += 4.0*sumA + 2.0*sumB;
	return sum*dx/3.0;
}
double CompositeSimpson(double** const &fEq, int const &End,int const &j,double const &dx)
{
	double sum = fEq[0][j] - fEq[End-1][j];
	double sumA = 0.0,sumB = 0.0;
	for(int i = 1;i < End-1;i+=2)
	{
		sumA += fEq[i][j];
		sumB += fEq[i+1][j];
	}
	sum += 4.0*sumA + 2.0*sumB;
	return sum*dx/3.0;
}
double CompositeSimpson(double* const &fEq,int const &End,double* const &xi,double const &dx)
{
	double sum = fEq[0]*xi[0] - fEq[End-1]*xi[End-1];
	double sumA = 0.0,sumB = 0.0;
	for(int i = 1;i < End-1;i+=2)
	{
		sumA += fEq[i]*xi[i];
		sumB += fEq[i+1]*xi[i+1];
	}
	sum += 4.0*sumA + 2.0*sumB;
	return sum*dx/3.0;
}
double CompositeSimpson(double* const &fEq,int const &End,
	double* const &xi_x,double* const &xi_y,double const &dx)
{
	double sum = fEq[0]*xi_x[0]*xi_y[0] - fEq[End-1]*xi_x[End-1]*xi_y[End-1];
	double sumA = 0.0,sumB = 0.0;
	for(int i = 1;i < End-1;i+=2)
	{
		sumA += fEq[i]*xi_x[i]*xi_y[i];
		sumB += fEq[i+1]*xi_x[i+1]*xi_y[i+1];
	}
	sum += 4.0*sumA + 2.0*sumB;
	return sum*dx/3.0;
}
//--------------------------------------Cotes-----------------------------------------------------
double CompositeCotes(double* const &fEq,int const &End,double const &dx)
{
	double sum = 3.5*(fEq[0] - fEq[End-1]);
	double sumA = 0.0,sumB = 0.0,sumC = 0.0;
	for(int i = 1;i < End-1;i+=4)
	{
		sumA += fEq[i] + fEq[i+2];
		sumB += fEq[i+1];
		sumC += fEq[i+3];
	}
	sum += 16.0*sumA + 6.0*sumB + 7.0*sumC;
	return sum*dx*4.0/45.0;
}
double CompositeCotes(double** const &fEq,int const &End,int const &j,double const &dx)
{
	double sum = 3.5*(fEq[0][j] - fEq[End-1][j]);
	double sumA = 0.0,sumB = 0.0,sumC = 0.0;
	for(int i = 1;i < End-1;i+=4)
	{
		sumA += fEq[i][j] + fEq[i+2][j];
		sumB += fEq[i+1][j];
		sumC += fEq[i+3][j];
	}
	sum += 16.0*sumA + 6.0*sumB + 7.0*sumC;
	return sum*dx*4.0/45.0;
}
double CompositeCotes(double* const &fEq,int const &End,double* const &xi,double const &dx)
{
	double sum = 3.5*(fEq[0]*xi[0] - fEq[End-1]*xi[End-1]);
	double sumA = 0.0,sumB = 0.0,sumC = 0.0;
	for(int i = 1;i < End-1;i+=4)
	{
		sumA += fEq[i]*xi[i] + fEq[i+2]*xi[i+2];
		sumB += fEq[i+1]*xi[i+1];
		sumC += fEq[i+3]*xi[i+3];
	}
	sum += 16.0*sumA + 6.0*sumB + 7.0*sumC;
	return sum*dx*4.0/45.0;
}
double CompositeCotes(double* const &fEq,int const &End,
	double* const &xi_x,double* const &xi_y,double const &dx) 
{
	double sum = 3.5*(fEq[0]*xi_x[0]*xi_y[0] - fEq[End-1]*xi_x[End-1]*xi_y[End-1]);
	double sumA = 0.0,sumB = 0.0,sumC = 0.0;
	for(int i = 1;i < End-1;i+=4)
	{
		sumA += fEq[i]*xi_x[i]*xi_y[i] + fEq[i+2]*xi_x[i+2]*xi_y[i+2];
		sumB += fEq[i+1]*xi_x[i+1]*xi_y[i+1];
		sumC += fEq[i+3]*xi_x[i+3]*xi_y[i+3];
	}
	sum += 16.0*sumA + 6.0*sumB + 7.0*sumC;
	return sum*dx*4.0/45.0;
}
#endif