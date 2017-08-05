#include "Matrix.h"

Matrix Matrix::operator*(Matrix mat1)
	{
	Matrix A(rows, mat1.getCols());
	if (mat1.getRows() != cols)
		{
		cout << "¾ØÕó²»ÄÜÏà³Ë" << endl;
		}
	else
		{
		for (int i = 0; i < rows; i++)
			{
			for (int j = 0; j <mat1.getCols(); j++)
				{
				for (int k = 0; k < cols; k++)
					{
					(A.mat)[i][j] += (mat[i][k] * (mat1.mat)[k][j]);
					}
				//cout << (A.mat)[i][j] << "  ";
				}
			}
		}
	return A;
	//return Matrix();
	}

Matrix Matrix::operator*(double c)
	{
	Matrix B(rows, cols);
	for (int i = 0; i < rows; i++)
		{
		for (int j = 0; j < cols; j++)
			{
			(B.mat)[i][j] = mat[i][j] * c;
			}
		}

	return B;
	}
Matrix Matrix::operator+(Matrix mat1)
	{
	Matrix T(rows, cols);
	for (int i = 0; i < rows; i++)
		{
		for (int j = 0; j < cols; j++)
			{
			(T.mat)[i][j] = this->GetMatrixElement(i, j) + (mat1.mat)[i][j];
			}
		}
	return T;
	}
Matrix Matrix::operator-(Matrix mat1)
	{
	Matrix T(rows, cols);
	for (int i = 0; i < rows; i++)
		{
		for (int j = 0; j < cols; j++)
			{
			(T.mat)[i][j] = this->GetMatrixElement(i, j) - (mat1.mat)[i][j];
			}
		}
	return T;
	}
Matrix Matrix::Transpose()
	{
	Matrix b(cols, rows);
	for (int i = 0; i < rows; i++)
		{
		for (int j = 0; j < cols; j++)
			{
			(b.mat)[j][i] = mat[i][j];
			}
		}
	return b;
	}
Matrix::~Matrix()
	{
	}
