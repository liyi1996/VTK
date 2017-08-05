#include"SparseMatrix.h"
SparseMatrix::SparseMatrix(const vector<vector<double>>& mat1)
	{
	int rows = mat1.size();//获取行数
	int cols = mat1[0].size();//获取列数
	this->rows = rows;
	this->cols = cols;
	int count1 = 0;//矩阵存储位置信息，比如第一个元素为0，第二个元素为1
	for (int i = 0; i < rows; i++)
		{
		rowsBegin.push_back(count1);
		for (int j = 0; j < cols; j++)
			{
			if (mat1[i][j] != 0)
				{
				matrixInformation.push_back(i);
				matrixInformation.push_back(j);
				matrixInformation.push_back(mat1[i][j]);
				rowsInColumns[i].push_back(j);
				columnsInRows[j].push_back(i);
				count1++;
				}
			}
		}
	}
SparseMatrix::SparseMatrix(const vector<double>&vec1)
	{
	int rows = 1;
	int cols = vec1.size();
	this->rows = rows;
	this->cols = cols;
	int count1 = 0;
	for (int i = 0; i < rows; i++)
		{
		rowsBegin.push_back(count1);
		for (int j = 0; j < cols; j++)
			{
			if (vec1[j] != 0)
				{
				matrixInformation.push_back(0);
				matrixInformation.push_back(j);
				matrixInformation.push_back(vec1[j]);
				rowsInColumns[0].push_back(j);
				columnsInRows[j].push_back(0);
				count1++;
				}
			}
		}
	}

SparseMatrix  SparseMatrix::operator+(const SparseMatrix &h)
	{
	SparseMatrix Result;
	int row = h.rows;
	int col = h.cols;
	Result.rows = row;
	Result.cols = col;
	vector<int>noZeroColumns;
	vector<int>noZeroColumnsInh;
	vector<int>::iterator it1, it2;
	int count1 = 0;
	for (int i = 0; i < row; i++)
		{
		noZeroColumnsInh = (h.rowsInColumns).at(i);
		noZeroColumns = getNoZeroColumns(i);
		it1 = noZeroColumnsInh.begin();
		it2 = noZeroColumnsInh.end();
		Result.rowsBegin.push_back(count1);
		//cout << noZeroColumns.size() << endl;
		for (int j = 0; j < noZeroColumns.size(); j++)
			{

				while (it1!=it2&&*it1 < noZeroColumns[j])
					{
					Result.matrixInformation.push_back(i);
					Result.matrixInformation.push_back(*it1);
					Result.matrixInformation.push_back(h.matrixInformation[3 * h.rowsBegin[i] + 3 * (it1 - noZeroColumnsInh.begin()) + 2]);
					Result.rowsInColumns[i].push_back(*it1);
					Result.columnsInRows[*it1].push_back(i);
					//cout << "*it1 samller" << endl;
					//cout << Result.matrixInformation.back() << endl;
					it1++;
					count1++;
					}
				if (it1 != it2&&*it1 == noZeroColumns[j])
					{
					Result.rowsInColumns[i].push_back(*it1);
					Result.columnsInRows[*it1].push_back(i);
					Result.matrixInformation.push_back(i);
					Result.matrixInformation.push_back(*it1);
					Result.matrixInformation.push_back(h.matrixInformation[3 * h.rowsBegin[i] + 3 * (it1 - noZeroColumnsInh.begin()) + 2] + matrixInformation[3 * rowsBegin[i] + 3 * j + 2]);
					it1++;
					count1++;
				//	cout << "*it1 equal  " << endl;
					//cout << Result.matrixInformation.back() << endl;
					}
				else if (it1!=it2&&*it1 > noZeroColumns[j])
					{
					Result.matrixInformation.push_back(i);
					Result.matrixInformation.push_back(noZeroColumns[j]);
					Result.matrixInformation.push_back(matrixInformation[3 * rowsBegin[i] + 3 * j + 2]);
					Result.rowsInColumns[i].push_back(noZeroColumns[j]);
					Result.columnsInRows[noZeroColumns[j]].push_back(i);
					count1++;
				//	cout << "*it1 bigger " << endl;
				//	cout << Result.matrixInformation.back() << endl;
					}
				else if (it1 == it2)
					{
					Result.matrixInformation.push_back(i);
					Result.matrixInformation.push_back(noZeroColumns[j]);
					Result.matrixInformation.push_back(matrixInformation[3 * rowsBegin[i] + 3 * j + 2]);
					Result.rowsInColumns[i].push_back(noZeroColumns[j]);
					Result.columnsInRows[noZeroColumns[j]].push_back(i);
					count1++;
					}
			
			}
		while (it1 != it2)
			{
			Result.matrixInformation.push_back(i);
			Result.matrixInformation.push_back(*it1);
			Result.matrixInformation.push_back(h.matrixInformation[3 * h.rowsBegin[i] + 3 * (it1 - noZeroColumnsInh.begin()) + 2]);
			Result.rowsInColumns[i].push_back(*it1);
			Result.columnsInRows[*it1].push_back(i);
			count1++;
			it1++;
			}
		}
	return Result;
	}

SparseMatrix SparseMatrix::operator*(const SparseMatrix &SP)
	{
	SparseMatrix result;
	result.rows = rows;
	result.cols = SP.cols;
	int count1 = 0;
	vector<int>noZeroColumns;
	vector<int>noZeroColumns2;
	vector<int>::iterator it;
	for (int i = 0; i < rows; i++)
		{
		result.rowsBegin.push_back(count1);
		for (int j = 0; j < SP.cols; j++)
			{
			
			double sum = 0;//每个元素之和
			noZeroColumns = getNoZeroColumns(i);

			for (int k = 0; k < noZeroColumns.size(); k++)
				{
				if (SP.rowsInColumns.find(noZeroColumns[k])!=SP.rowsInColumns.end())
					noZeroColumns2 = SP.rowsInColumns.at(noZeroColumns[k]);
				else
					continue;
				//cout << noZeroColumns2.size() << endl;
				if (!noZeroColumns2.empty())
					{
					if (find(noZeroColumns2.begin(), noZeroColumns2.end(), j) != noZeroColumns2.end())
						{
						it = find(noZeroColumns2.begin(), noZeroColumns2.end(), j);
						sum += matrixInformation[3 * rowsBegin[i] + 3 * k + 2] * SP.matrixInformation[3 * SP.rowsBegin[noZeroColumns[k]] + 3 * (it - noZeroColumns2.begin()) + 2];
						}
					}
				}
			if (sum != 0)
				{
				result.matrixInformation.push_back(i);
				result.matrixInformation.push_back(j);
				result.matrixInformation.push_back(sum);
				result.columnsInRows[j].push_back(i);
				result.rowsInColumns[i].push_back(j);
				count1++;
				}
			}
		}
	
	return result;
	}

SparseMatrix SparseMatrix::operator-(const SparseMatrix &A)
	{
	SparseMatrix result;
	SparseMatrix Atem = A;
	for (int i = 2; i < Atem.matrixInformation.size(); i += 3)
		{
		Atem.matrixInformation[i] = -A.matrixInformation[i];
		}
	result = *this + Atem;
	return result;
	}

SparseMatrix SparseMatrix::operator*(double c)
	{
	SparseMatrix result;
	result = *this;
	for (int i = 2; i <result.matrixInformation.size(); i += 3)
		{
		//cout << "i:  " << i << endl;
		//cout << matrixInformation.size() << endl;
		result.matrixInformation[i] =c*result.matrixInformation[i];
		}
	return result;
	}

SparseMatrix SparseMatrix::Transpose()
	{
	SparseMatrix Result;
	int count1 = 0;
	Result.cols = rows;
	Result.rows = cols;
	vector<int> noZeroRows;
	vector<int> noZeroColumns;
	vector<int>::iterator it;
	for (int j = 0; j < cols; j++)
		{
		Result.rowsBegin.push_back(count1);
		noZeroRows = getNoZeroRows(j);
		for (int i = 0; i < noZeroRows.size(); i++)
			{
			Result.matrixInformation.push_back(j);
			Result.matrixInformation.push_back(noZeroRows[i]);
		//cout << noZeroRows[i] << endl;
		//cout << rowsBegin[2] << endl;
		noZeroColumns = getNoZeroColumns(noZeroRows[i]);
		it = find(noZeroColumns.begin(), noZeroColumns.end(), j);

			Result.matrixInformation.push_back(matrixInformation[3 * rowsBegin[noZeroRows[i]] + 3 * (it-noZeroColumns.begin()) + 2]);
			Result.rowsInColumns[j].push_back(noZeroRows[i]);
			Result.columnsInRows[noZeroRows[i]].push_back(j);
			count1++;
			}
		}
	return Result;
	}
