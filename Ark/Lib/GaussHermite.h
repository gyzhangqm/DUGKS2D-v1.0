double IntegralGH(double* const &f,int Num,double * const &wGH)
{
	double sum = 0;
	for(int i = 0;i < Num;++i)
	{
		sum += f[i]*wGH[i];
	}
	return sum;
}
double IntegralGH(double** const &f,int Num,int const &j,double * const &wGH)
{
	double sum = 0;
	for(int i = 0;i < Num;++i)
	{
		sum += f[i][j]*wGH[i];
	}
	return sum;
}
double IntegralGH(double * const &f,int Num,double * const &xiA,double * const &wGH)
{
	double sum = 0;
	for(int i = 0;i < Num;++i)
	{
		sum += xiA[i]*f[i]*wGH[i];
	}
	return sum;
}
double IntegralGH(double * const &f,int Num,
					double * const &xiA,double * const &xiB,double * const &wGH)
{
	double sum = 0;
	for(int i = 0;i < Num;++i)
	{
		sum += xiA[i]*xiB[i]*f[i]*wGH[i];
	}
	return sum;
}