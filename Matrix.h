#pragma once
#include<iostream>
#include<vector>
using namespace std;
class Matrix
	{

	public:
		Matrix()
			{

			}
		Matrix(int row, int col) :rows(row), cols(col) {
			mat.resize(row);
			for (int i = 0; i < row; i++)
				{
				mat[i].resize(col);
				}
			};
		Matrix(vector<double> ans)
			{
			rows = 1;
			cols = ans.size();
			//cout << "cols " << cols << endl;
			mat.resize(1);
			mat[0].resize(cols);
			for (int i = 0; i < cols; i++)
				{
				mat[0][i] = ans[i];
				}
			}
		Matrix(const Matrix&) = default;
		void PrintMatrix() const
			{
			for (int i = 0; i < rows; i++)
				{
				for (int j = 0; j < cols; j++)
					{
					cout << this->GetMatrixElement(i, j) << "  ";
					}
				cout << endl;
				}
			}
		void SetMatrix(vector<vector<double> > Mat)//行列必须相等
			{
			mat = Mat;
			rows = Mat.size();
			cols = Mat[0].size();
			}
		void SetMatrixElement(int row, int col, double r)
			{
			mat[row][col] = r;
			}
		double GetMatrixElement(int row, int col) const
			{
			return mat[row][col];
			}
		void setRows(int row)
			{
			rows = row;
			}
		int getRows() const
			{
			return rows;
			}
		void setCols(int col)
			{
			cols = col;
			}
		int getCols()const
			{
			return cols;
			}
		Matrix operator*(Matrix);
		Matrix operator*(double);
	//	Matrix operator*(double)
		vector<double> getRowVectors() const
			{
			if (rows == 1)
			return mat[0];
			}
		Matrix operator+(Matrix);
		Matrix operator-(Matrix);
		Matrix Transpose();
		double VectorTMultipleMatrixVector(Matrix & mat2)
			{
			//mat2.PrintMatrix();
			if (mat2.getRows() == 1)
				{
				Matrix middle(1, mat2.getCols());
				//cout << cols << endl;
				for (int i = 0; i < cols; i++)
					{
					double k=0;
					for (int j = 0; j < cols; j++)
						{
						k += mat2.GetMatrixElement(0, j)* this->GetMatrixElement(j, i);
						//cout << this->GetMatrixElement(j, i) << endl;
						}
					middle.SetMatrixElement(0, i, k);
			
					}
				//middle.PrintMatrix();
				//PrintMatrix();
				double sum = 0;
				for (int i = 0; i < cols; i++)
					{
					sum += middle.GetMatrixElement(0, i)*mat2.GetMatrixElement(0, i);
					}
				
				return sum;
				}
			else if (mat2.getCols() == 1)
				{
				Matrix middle(1, mat2.getRows());
				for (int i = 0; i < rows; i++)
					{
					double k = 0;
					for (int j = 0; j < rows; j++)
						{
						k += mat2.GetMatrixElement(j, 0)*this->GetMatrixElement(j, i);
						}
					middle.SetMatrixElement(0, i, k);
					}
				double sum = 0;
				for (int i = 0; i < rows; i++)
					{
					sum += middle.GetMatrixElement(0, i)*mat2.GetMatrixElement(i,0);
					}
				return sum;
				}
			}
		~Matrix();
	private:
		int rows=0;
		int cols=0;
		vector<vector<double> > mat;
	};