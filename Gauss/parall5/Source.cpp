// Gauss.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>

void PrintMatrix(double **matrix, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size + 1; j++)
		{

			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Gaus(double **matrix, int size)
{
	int row = 0, max_el_row;
	double max_el = 0.0, zero1 = 0.0, zero2 = 0.0, zero3 = 0.0;

	for (int j = 0; j < size; j++)  //cтолб
	{
#pragma omp parallel for
		for (int h = j; h < size; h++) //строка
		{
			if (fabs(matrix[h][j]) > max_el)
			{
				max_el_row = h;
				max_el = fabs(matrix[h][j]);
			}
		}

		max_el = 0.0;

		std::swap(matrix[max_el_row], matrix[row]);

		PrintMatrix(matrix, size);

#pragma omp parallel for
		for (int k = 1 + row; k < size; k++)  
		{
			zero1 = matrix[k][j] / matrix[row][j];

			for (int i = 0; i <= size; i++)
			{
				if (fabs(zero1) == 0)
					break;

				zero2 = (matrix[row][i] * zero1);
				zero3 = matrix[k][i] - zero2;
				matrix[k][i] = zero3;
			}
		}
		row++;
	}
}

void GausX(double **matrix, double *x, int size)
{
#pragma omp parallel for
	for (int i = size - 1; i >= 0; i--)
	{
		double s = 0;
		
		for (int j = size - 1; j > i; j--)
		{
			s += x[j] * matrix[i][j];
		}

		x[i] = (matrix[i][size] - s) / matrix[i][i];

	}
}

int main()
{
	srand(time(0));
	int size = 6;
	double **matrix = new double*[size];
	for (int i = 0; i < size; i++)
		matrix[i] = new double[size + 1];

	for (auto i = 0; i < size; i++)
		for (auto j = 0; j < size + 1; j++)
		{
			matrix[i][j] = rand() % 10 + 1;
		}

	Gaus(matrix, size);

	PrintMatrix(matrix, size);

	double *x = new double[size];

	GausX(matrix, x, size);

	for (int i = 0; i < size; i++)
		std::cout << x[i] << " ";
	std::cout << std::endl;

	system("pause");
	return 0;
}

