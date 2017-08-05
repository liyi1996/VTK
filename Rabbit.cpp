#include"Rabbit.h"
#include"Matrix.h"
void Rabbit:: CalcForce()
	{
	force.resize(2 * points);//force[0]代表Fx1 force[1]代表Fy1
	force.assign(2 * points, 0);
	for (int i = 0; i < 2 * points; i++)
		{
		if (Constraint[i / 2] == 1)
			{
			force[i] = 0;
			continue;
			}
		if (i % 2 == 0)
			{
			force[i] -= m*9.8;//重力
			}
		for (int j = 0; j < vtkevc[i / 2]->GetNumberOfIds(); j++)
			{
			int ConnectPoint = vtkevc[i / 2]->GetId(j);
			if (i % 2 == 0)//x坐标
				{
				double distance = DistanceBetwweenTwoPoints(i / 2, ConnectPoint);
				force[i] += stiff*(getPointX(ConnectPoint) - getPointX(i / 2)) / distance *(distance - length[i / 2][ConnectPoint]);
				}
			else
				{
				double distance = DistanceBetwweenTwoPoints(i / 2, ConnectPoint);
				force[i] += stiff*(getPointY(ConnectPoint) - getPointY(i / 2)) / distance *(distance - length[i / 2][ConnectPoint]);
				}


			}
		}
	}
void Rabbit::CalcNextForce()
	{
	forcenext.resize(2 * points);//force[0]代表Fx1 force[1]代表Fy1
	forcenext.assign(2 * points, 0);
	for (int i = 0; i < 2 * points; i++)
		{
		if (Constraint[i / 2] == 1)
			continue;
		if (i % 2 == 0)
			{
			forcenext[i] -= m*9.8;//重力
			}
		for (int j = 0; j < vtkevc[i / 2]->GetNumberOfIds(); j++)
			{
			int ConnectPoint = vtkevc[i / 2]->GetId(j);
			if (i % 2 == 0)//x坐标
				{
				double distance = sqrt((locationnext[2 * ConnectPoint] - locationnext[i])*(locationnext[2 * ConnectPoint] - locationnext[i]) + (locationnext[2 * ConnectPoint + 1] - locationnext[i + 1])*(locationnext[2 * ConnectPoint + 1] - locationnext[i + 1]));
				forcenext[i] += stiff*(locationnext[2 * ConnectPoint] - locationnext[i]) / distance *(distance - length[i / 2][ConnectPoint]);
				}
			else
				{
				double distance = sqrt((locationnext[2 * ConnectPoint] - locationnext[i - 1])*(locationnext[2 * ConnectPoint] - locationnext[i - 1]) + (locationnext[2 * ConnectPoint + 1] - locationnext[i])*(locationnext[2 * ConnectPoint + 1] - locationnext[i]));
				forcenext[i] += stiff*(locationnext[2 * ConnectPoint + 1] - locationnext[i]) / distance *(distance - length[i / 2][ConnectPoint]);
				}
			}
		}
	}
vtkSmartPointer<vtkIdList> Rabbit:: GetConnectedVertices(int id)
	{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();

	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	this->rabbitInformation->GetPointCells(id, cellIdList);

	/*
	cout << "Vertex 0 is used in cells ";
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
	cout << cellIdList->GetId(i) << ", ";
	}
	cout << endl;
	*/

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
		{
		//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		this->rabbitInformation->GetCellPoints(cellIdList->GetId(i), pointIdList);
		//cout << "pointIdList  " << endl;
		//cout << pointIdList->GetNumberOfIds() << endl;
		//cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;

		if (pointIdList->GetId(0) != id)
			{
			connectedVertices->InsertUniqueId(pointIdList->GetId(0));
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			}
		if (pointIdList->GetId(1) != id)
			{
			connectedVertices->InsertUniqueId(pointIdList->GetId(1));
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			}
		if (pointIdList->GetId(2) != id)
			{
			connectedVertices->InsertUniqueId(pointIdList->GetId(2));
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			}
		}

	return connectedVertices;
	}
vector<double> ConjugateGradientMethod(vector<vector<double> >A, vector<double> B)
	{
	vector<double> ans;
	ans.resize(B.size());
	Matrix AA;
	AA.SetMatrix(A);
	//SparseMatrix AA(A);
	Matrix BB(B);
	//SparseMatrix BB(B);
	vector<double> x0;
	x0.resize(B.size());
	Matrix X0(x0);
	//SparseMatrix X0(x0);
	//cout << X0.getRows() << endl;
	vector<double> x1;
	x1.resize(B.size());
	Matrix X1(x1);
	//SparseMatrix X1(x1);
	vector<double> r0;
	r0.resize(B.size());
	Matrix R0(r0);
	//SparseMatrix R0(r0);
	vector<double> p0;
	p0.resize(B.size());
	Matrix P0(p0);
	//SparseMatrix P0(p0);
	vector<double> p1;
	p1.resize(B.size());
	Matrix P1(p1);
	//SparseMatrix P1(p1);
	vector<double> r1;
	r1.resize(B.size());
	Matrix R1(r1);
	//SparseMatrix R1(r1);
	double alpha;
	double beta;
	//X0.Transpose().printInformation();
	//BB.printInformation();
	//if (!X0.empty())
		R0 = BB - (AA*X0.Transpose()).Transpose();
	//else
		//R0 = BB;
	P0 = R0;
	//Matrix PP(p0);
	double t = 0;
	double tdown = 0;
	//SparseMatrix SS;
	while (1)
		{
		//计算rk的转置*rk
		tdown = 0;
		t = 0;
		t = (R0*R0.Transpose()).GetMatrixElement(0, 0);
		//R0.printInformation();
		//R0.Transpose().printInformation();
		//计算pk的转置乘以A矩阵乘以pk
	//	cout << "t" << t << endl;
	     
		//SS.printInformation();
	
		tdown = (P0*AA*P0.Transpose()).GetMatrixElement(0,0);
		//cout << tdown << endl;
		alpha = t / tdown;
		//计算xk+1
		X1 = X0 + P0*alpha;
		//X1.printInformation();
		R1 = R0 - (AA*P0.Transpose()*alpha).Transpose();
		double exit = 0.0;//判断结束条件，计算rk+1的模
		exit = (R1*R1.Transpose()).GetMatrixElement(0, 0);
	//	cout << exit << endl;
		if (exit < 10E-3)
			{
			ans = X1.getRowVectors();
			return ans;
			}
		//计算betak
		beta = exit / t;
		P1 = R1 + P0*beta;
		X0 = X1;
		P0 = P1;
		R0 = R1;
		}
	}
void Rabbit::CalcpartialF()
	{
	partialF.resize(2 * points);
	for (int i = 0; i < 2 * points; i++)
		{
		partialF[i].resize(2 * points);//初始化偏导数数组
		partialF[i].assign(2 * points, 0);
		}

	for (int i = 0; i < 2 * points; i++)
		{
		if (Constraint[i / 2] == 1)
			{
			continue;
			}
		for (int j = 0; j < vtkevc[i / 2]->GetNumberOfIds(); j++)
			{
			int ConnectPoint = vtkevc[i / 2]->GetId(j);

			if (i % 2 == 0)//x坐标
				{
				double distance = DistanceBetwweenTwoPoints(i / 2, ConnectPoint);
				double Fxx = -stiff*(distance - length[i / 2][ConnectPoint]) / distance - stiff*(getPointX(ConnectPoint) - getPointX(i / 2))*(getPointX(ConnectPoint) - getPointX(i / 2)) / distance / distance / distance
					*(length[i / 2][ConnectPoint]);
				double Fxy = stiff*(getPointX(ConnectPoint) - getPointX(i / 2))*(-(getPointY(ConnectPoint) - getPointY(i / 2))) / distance / distance / distance *(length[i / 2][ConnectPoint]);
				//partialF[i][i] += -stiff*(distance - length[i / 2][ConnectPoint]) / distance - stiff*(getPointX(ConnectPoint) - getPointX(i / 2))*(getPointX(ConnectPoint) - getPointX(i / 2)) / distance / distance / distance
				//*(length[i / 2][ConnectPoint]);
				//partialF[i][i + 1] += stiff*(getPointX(ConnectPoint) - getPointX(i / 2))*(-(getPointY(ConnectPoint) - getPointY(i / 2))) / distance / distance / distance *(length[i / 2][ConnectPoint]);
				//partialF[i][2 * ConnectPoint] += (-(-stiff*(distance - length[i / 2][ConnectPoint]) / distance  - stiff*(getPointX(ConnectPoint) - getPointX(i / 2))*(getPointX(ConnectPoint) - getPointX(i / 2)) / distance / distance / distance 
				//*(length[i / 2][ConnectPoint])));
				//partialF[i][2 * ConnectPoint+1] += (-(stiff*(getPointX(ConnectPoint) - getPointX(i / 2))*(-(getPointY(ConnectPoint) - getPointY(i / 2))) / distance / distance / distance*( length[i / 2][ConnectPoint])));
				partialF[i][i] += Fxx;
				partialF[i][i + 1] += Fxy;
				partialF[i][2 * ConnectPoint] += (-Fxx);
				partialF[i][2 * ConnectPoint + 1] += (-Fxy);
				}
			else
				{
				double distance = DistanceBetwweenTwoPoints(i / 2, ConnectPoint);
				double Fyy = -stiff*(distance - length[i / 2][ConnectPoint]) / distance - stiff*(getPointY(ConnectPoint) - getPointY(i / 2))*(getPointY(ConnectPoint) - getPointY(i / 2)) / distance / distance / distance
					*(length[i / 2][ConnectPoint]);
				double Fyx = stiff*(getPointY(ConnectPoint) - getPointY(i / 2))*(-(getPointX(ConnectPoint) - getPointX(i / 2))) / distance / distance / distance*(length[i / 2][ConnectPoint]);

				//	partialF[i][i] += -stiff*(distance - length[i / 2][ConnectPoint]) / distance  - stiff*(getPointY(ConnectPoint) - getPointY(i / 2))*(getPointY(ConnectPoint) - getPointY(i / 2)) / distance / distance / distance
				//	*( length[i / 2][ConnectPoint] );
				//partialF[i][i -1] += stiff*(getPointY(ConnectPoint) - getPointY(i / 2))*(-(getPointX(ConnectPoint) - getPointX(i / 2))) / distance / distance / distance*(length[i / 2][ConnectPoint]);
				//partialF[i][2 * ConnectPoint] += (-(stiff*(getPointY(ConnectPoint) - getPointY(i / 2))*(-(getPointX(ConnectPoint) - getPointX(i / 2))) / distance / distance / distance *( length[i / 2][ConnectPoint])));
				//partialF[i][2 * ConnectPoint+1] += (-(-stiff*(distance - length[i / 2][ConnectPoint]) / distance  - stiff*(getPointY(ConnectPoint) - getPointY(i / 2))*(getPointY(ConnectPoint) - getPointY(i / 2)) / distance / distance / distance 
				//*(length[i / 2][ConnectPoint] )));
				partialF[i][i] += Fyy;
				partialF[i][i - 1] += Fyx;
				partialF[i][2 * ConnectPoint] += (-Fyx);
				partialF[i][2 * ConnectPoint + 1] += (-Fyy);
				}
			}
		}
	}