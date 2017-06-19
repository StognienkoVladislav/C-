#include <math.h>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>

#define A_DIFF 1.0
#define B_DIFF 1.0
#define C_DIFF 1.0
#define A_FUNC 1.0
#define B_FUNC 1.0

#define X_MIN 0.0
#define X_MAX 1.0
#define X_POINT_AMOUNT 10
#define X_STEP ((X_MAX - X_MIN)/(X_POINT_AMOUNT - 1))

#define T_MIN 0.0
#define T_MAX 1.0
#define T_POINT_AMOUNT 100
#define T_STEP ((T_MAX - T_MIN)/(T_POINT_AMOUNT - 1))

long double xtStepFunction(long double wiPlus1, long double wiMinus1, long double wi, long double tau, long double h)
{
	return wi + tau * (A_DIFF * (log(wi)*((wiPlus1 - wiMinus1) / (2 * h))*((wiPlus1 - wiMinus1) / (2 * h))
		+ (1.0 / wi) * ((wiMinus1 - 2 * wi + wiPlus1) / (h * h)))
		+ B_DIFF + C_DIFF * wi);
}
long double xtCorrectFunction(long double x, long double t)
{
	return A_FUNC * exp(C_DIFF * t - (B_DIFF * x * x) / (2 * A_DIFF) + B_FUNC * x);
}
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int myid, numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); //число процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//номер процесса

	long double matrixOfApproximateAnswers[T_POINT_AMOUNT][X_POINT_AMOUNT] = { 0 };
	long double matrixOfCorrectAnswers[T_POINT_AMOUNT][X_POINT_AMOUNT] = { 0 };

	if (myid == 0) //главн процесс  
	{
#pragma omp parallel for //матрица правильных значений 
		for (int i = 0; i < T_POINT_AMOUNT; i++) {
			for (int j = 0; j < X_POINT_AMOUNT; j++) {
				matrixOfCorrectAnswers[i][j] = xtCorrectFunction(X_MIN + j * X_STEP, T_MIN + i * T_STEP);
			}
		}
		// Algorithm starts 
#pragma omp parallel for //заполнить 0 шар
		for (int i = 0; i < X_POINT_AMOUNT; i++)
		{
			matrixOfApproximateAnswers[0][i] = xtCorrectFunction(X_MIN + i * X_STEP, T_MIN);
		}
#pragma omp parallel for // хmin и хmax
		for (int i = 0; i < T_POINT_AMOUNT; i++)
		{
			matrixOfApproximateAnswers[i][0] = xtCorrectFunction(X_MIN, T_MIN + i * T_STEP);
			matrixOfApproximateAnswers[i][X_POINT_AMOUNT - 1] = xtCorrectFunction(X_MAX, T_MIN + i * T_STEP);
		}
	}
	for (int i = 1; i < T_POINT_AMOUNT; i++)// заполнить 1 шар 
	{
		MPI_Bcast(matrixOfApproximateAnswers[i - 1], X_POINT_AMOUNT, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
#pragma omp parallel for
		for (int j = myid + 1; j < X_POINT_AMOUNT - 1; j += numprocs)
		{
			matrixOfApproximateAnswers[i][j] = xtStepFunction(matrixOfApproximateAnswers[i - 1][j + 1],
				matrixOfApproximateAnswers[i - 1][j - 1], matrixOfApproximateAnswers[i - 1][j], T_STEP, X_STEP);
		}
		if (myid == 0)// собирает все данн
		{
			for (int k = 1; k < numprocs; k++)
			{
				long double temp[X_POINT_AMOUNT];
				MPI_Recv(temp, X_POINT_AMOUNT, MPI_LONG_DOUBLE, k, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int j = k; j < X_POINT_AMOUNT - 1; j += numprocs)
				{
					matrixOfApproximateAnswers[i][j] = temp[j];
				}
			}
		}
		else//отправл€ют 
		{
			MPI_Send(&matrixOfApproximateAnswers[i], X_POINT_AMOUNT, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
	if (myid == 0)
	{
		// Saving results to file
		std::ofstream fout1("results_clear.dat");
		std::ofstream fout2("results_omp.dat");
		for (int i = 0; i < T_POINT_AMOUNT; i++)
		{
			for (int j = 0; j < X_POINT_AMOUNT; j++)
			{
				fout1 << matrixOfCorrectAnswers[i][j] << '\t';
				fout2 << matrixOfApproximateAnswers[i][j] << '\t';
			}
			fout1 << std::endl;
			fout2 << std::endl;
		}
		fout1.close();
		fout2.close();
	}

	MPI_Finalize();
	return 0;
}
