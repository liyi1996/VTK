#include"Matrix.h"
#include"SparseMatrix.h"
#include <vtkOBJReader.h>
#include<vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include<vector>
#include<cmath>
using namespace std;
vector<double> ConjugateGradientMethod(vector<vector<double> >A, vector<double> B);
class Rabbit
	{

	public:
		vtkSmartPointer<vtkPolyData> getPolydata() const
			{
			return rabbitInformation;
			}
		Rabbit(vtkSmartPointer<vtkPolyData> copy)
			{
			rabbitInformation = copy;
			points = rabbitInformation->GetNumberOfPoints();
			cells = rabbitInformation->GetNumberOfCells();
			m = 1;
			h =0.01;
			stiff = 10000;
			}
		vtkSmartPointer<vtkIdList> GetConnectedVertices(int id);
		void GetAllConnectVertices() 
			{
			for (int i = 0; i < points; i++)
				{
				vtkSmartPointer<vtkIdList> tem;
				tem = GetConnectedVertices(i);
				vtkevc.push_back(tem);
				}
			}
		void PrintAllConnectVertices() const
			{
			for (int i = 0; i < points; i++)
				{
				cout << "i:  " << i << endl;
				for (int j = 0; j < vtkevc[i]->GetNumberOfIds(); j++)
					{
					cout << vtkevc[i]->GetId(j) << "  ";
					}
				cout << endl;
				}
			}
		bool isConnectedVertices(int a, int b)const
			{
			for (int i = 0; i < vtkevc[a]->GetNumberOfIds(); i++)
				{
				if (vtkevc[a]->GetId(i) == b)
					return true;
				}
			return false;
			}
		double DistanceBetwweenTwoPoints(int a, int b) const
			{
			double ax1, ay1;
			double bx1, by1;
			ax1 = (rabbitInformation->GetPoint(a))[0];
			ay1= (rabbitInformation->GetPoint(a))[1];
			bx1 = (rabbitInformation->GetPoint(b))[0];
			by1 = (rabbitInformation->GetPoint(b))[1];
			return sqrt((ax1 - bx1)*(ax1 - bx1) + (ay1 - by1)*(ay1 - by1));
			}
		void SetLength()
			{
			length.resize(points);
			for (int i = 0; i < points; i++)
				{
				vector<double> temlength;
				for (int j = 0; j < points; j++)
					{
					if (isConnectedVertices(i, j))
						{
						//cout << "i:  " << i << "j:  " << j << endl;
						temlength.push_back(DistanceBetwweenTwoPoints(i, j));
						}
					else
						{
						temlength.push_back(0);
						}
					}
				//cout << endl;
				length[i] = temlength;
				//cout << length[i].size() << endl;
				temlength.clear();
				}
			}
		void PrintLength() const
			{
			for (int i = 0; i < points; i++)
				{
				for (int j = 0; j < points; j++)
					{
					cout << length[i][j] << "  ";
					}
				cout << endl;
				}
			}
		double getPointX(int id) const
			{
			return (rabbitInformation->GetPoint(id))[0];
			}
		double getPointY(int id) const
			{
			return (rabbitInformation->GetPoint(id))[1];
			}
		double getPointZ(int id) const
			{
			return (rabbitInformation->GetPoint(id))[2];
			}
		void CalcForce();
		void PrintForce() const
			{
			for (int i = 0; i < 2 * points; i++)
				cout << force[i] << "  ";
			cout << endl;
			}
		void CalcpartialF();
		void PrintpartialF() const
			{
			int count1 = 0;
			for (int i = 0; i < 2 * points; i++)
				{
				for (int j = 0; j < 2 * points; j++)
					{
					cout << partialF[i][j] << "  ";
					if (partialF[i][j] != 0)
						count1++;
					}
				cout << endl;
				}
			cout << "count:  " << count1 << endl;
			}
		void InitLocation()
			{
			location.resize(points * 2);
			for (int i = 0; i < points; i++)
				{
				location[2 * i] = getPointX(i);
				location[2 * i + 1] = getPointY(i);
				}
			}
		void InitSpeed()
			{
			speed.resize(2 * points);
			for (int i = 0; i < points; i++)
				{
				speed[2 * i] =0;
				speed[2 * i + 1] =0;
				}
			}
		vector<double> CalcMinusF() const
			{
			vector<double> ans;
			ans.resize(2 * points);
			ans.assign(2 * points, 0.0);
			for (int i = 0; i < 2 * points; i++)
				{
				ans[i] = (h*h / m*force[i]) +h*speed[i];
				}
			return ans;
			}
		vector<vector<double> >CalcEquationPartial() const
			{
			vector<vector<double> >ans;
			ans.resize(2 * points);
			for (int i = 0; i < 2 * points; i++)
				{
				ans[i].resize(2 * points);
				ans[i].assign(2 * points, 0);
				}
			for (int i = 0; i < 2 * points; i++)
				{
				for (int j = 0; j < 2 * points; j++)
					{
					ans[i][j] += -h*h / m*partialF[i][j];
					}
				ans[i][i] += 1;
				}
			return ans;
			}
		void CalcNextLocation()
				{
			locationnext.resize(2 * points);
			locationnext.assign(2 * points, 0);
			vector<double> c;
		   c = ConjugateGradientMethod(CalcEquationPartial(), CalcMinusF());	
			for (int i = 0; i < points * 2; i++)
				{
				if (Constraint[i / 2] == 1)
					{
					locationnext[i] = location[i];
					}
				else
					locationnext[i] = c[i] +1* location[i];
				
			
				}
			}
		void CalcNextForce();
		void CalcNextSpeed()
			{
			speednext.resize(2 * points);
			speednext.assign(2 * points, 0);
			for (int i = 0; i < points * 2; i++)
				{
				speednext[i] = speed[i] + h / m*forcenext[i];
				}
			//cout << endl;
			}
		void InitAgain()
			{
			location = locationnext;
			speed = speednext;
			for (int i = 0; i < rabbitInformation->GetNumberOfPoints(); i++)
				{
				double z = this->getPointZ(i);
				rabbitInformation->GetPoints()->SetPoint(i, location[2 * i], location[2 * i + 1], z);
				}
			}
		void SetConstraint()
			{
			Constraint.resize(points);
			Constraint.assign(points, 0);
			for (int i = 0; i < points; i++)
				{
				if (getPointY(i)<-8)
					{
					Constraint[i] = 1;
					}
				}
		
			}
		~Rabbit()
			{
			}
	private:
		vtkSmartPointer<vtkPolyData> rabbitInformation;//兔子信息
	    int points;//点数
		int cells;//三角形数
		double m;//每个点的质量
		double h;//每个点的步长
		double stiff;//劲度系数
		vector<vector<double> > length;//弹簧原长数组
		vector<double> force;//弹力数组（包括重力）
		vector<double> forcenext;//下一弹力数组（包括重力）
		vector<vector<double> >partialF;//偏导数数组；
		vector<vtkSmartPointer<vtkIdList> >vtkevc;//计算连接的顶点数
		vector<double> speed;//速度数组
		vector<double> location;//位置数组
		vector<double> speednext;//下一速度数组
		vector<double> locationnext;//下一位置数组
		vector<int> Constraint;//限制数组
	};