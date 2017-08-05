#pragma once
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<algorithm>
using namespace std;
class SparseMatrix
	{
	public:
		SparseMatrix() = default;//默认构造函数
		SparseMatrix(const vector<vector<double>>&);//初始化矩阵
		SparseMatrix(const vector<double>&);//初始化向量
		SparseMatrix operator+(const SparseMatrix&);
		SparseMatrix operator*(const SparseMatrix&);
		SparseMatrix operator-(const SparseMatrix&);
		SparseMatrix operator*(double);
		bool empty() const
			{
			if (matrixInformation.empty())
				return true;
			return false;
			}
		int getCols() const
			{
			return cols;
			}
		int getRows() const
			{
			return rows;
			}
		double getElement(int row, int col)
			{
			vector<int> nozeroColumns = rowsInColumns[row];
			vector<int>::iterator it;
			if (find(nozeroColumns.begin(), nozeroColumns.end(), col) != nozeroColumns.end())
				{
				it = find(nozeroColumns.begin(), nozeroColumns.end(), col);
				return matrixInformation[3 * rowsBegin[row] + 3 * (it - nozeroColumns.begin()) + 2];
				}
			else
				return 0.0;
			}
		vector<int> getNoZeroColumns(int row)
			{
			map<int, vector<int>>::iterator it = rowsInColumns.find(row);
			if (it != rowsInColumns.end())
				{
				return it->second;
				}
			else
				return vector<int>();
			}
		vector<int> getNoZeroRows(int col)
			{
			map<int, vector<int>>::iterator it = columnsInRows.find(col);
			if (it != columnsInRows.end())
				{
				return it->second;
				}
			else
				return vector<int>();
			}
		SparseMatrix Transpose();
		void printInformation() const
			{
			for (int i = 0; i < matrixInformation.size(); i++)
				{
				if (i % 3 == 0)
					cout << endl;
				cout << matrixInformation[i] << "  ";
				
				}
			cout << endl;
			}
		~SparseMatrix()=default;//默认析构函数
	private:
		vector<double> matrixInformation;//稀疏矩阵的存储信息 0代表行数 1代表列数 2代表元素值
		map<int, vector<int> >rowsInColumns;//每一行数对应的非零列数
		map<int, vector<int>>columnsInRows;//每一列数对应的非零行数
		int rows=0;//行数
		int cols=0;//列数
		vector<int>rowsBegin;//代表每一行的起始位置
	};