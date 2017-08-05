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
		SparseMatrix() = default;//Ĭ�Ϲ��캯��
		SparseMatrix(const vector<vector<double>>&);//��ʼ������
		SparseMatrix(const vector<double>&);//��ʼ������
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
		~SparseMatrix()=default;//Ĭ����������
	private:
		vector<double> matrixInformation;//ϡ�����Ĵ洢��Ϣ 0�������� 1�������� 2����Ԫ��ֵ
		map<int, vector<int> >rowsInColumns;//ÿһ������Ӧ�ķ�������
		map<int, vector<int>>columnsInRows;//ÿһ������Ӧ�ķ�������
		int rows=0;//����
		int cols=0;//����
		vector<int>rowsBegin;//����ÿһ�е���ʼλ��
	};