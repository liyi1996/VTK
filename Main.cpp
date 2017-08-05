#include"Matrix.h"
#include"SparseMatrix.h"
#include"Rabbit.h"
#include <vtkOBJReader.h>
#include<vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<Windows.h>
#include<ctime>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
using namespace std;
int main()
	{
	DWORD t1, t2;
	/*
	SparseMatrix A{ vector<double>{1,2,3} };
	SparseMatrix B = A.Transpose();
	SparseMatrix C(vector<vector<double>>{ vector<double>{1,0,4},vector<double>{2,1,5},vector<double>{3,2,6} });
	SparseMatrix D = A*C*B;
	D.printInformation();
	SparseMatrix E = A*C;
	E.printInformation();
	SparseMatrix F = E*B;
	F.printInformation();
	SparseMatrix G = F*9.8;
	G.printInformation();
	SparseMatrix H = G - F;
	H.printInformation();
	*/
	//cout << B.getCols() << "  " << B.getRows() << endl;
	//SparseMatrix C = A*B;
	//cout << C.getCols() << "  " << C.getRows() << endl;
	//C.printInformation();
	t1 = GetTickCount();
	string filename = "D:\\Download\\rabbit\\funny.obj";
	vtkSmartPointer<vtkOBJReader> reader =
		vtkSmartPointer<vtkOBJReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata = reader->GetOutput();

	
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polydata);
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetRepresentationToWireframe();
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green
	renderer->SetViewport(0, 0.1, 1, 0.9);
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(1900, 600);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	
	//renderWindow->Render();
	Rabbit rabbit(polydata);
	rabbit.GetAllConnectVertices();
	//rabbit.PrintAllConnectVertices();
	rabbit.SetConstraint();
	rabbit.SetLength();
//rabbit.PrintLength();
	rabbit.CalcForce();
	//rabbit.PrintForce();
	rabbit.CalcpartialF();
	//rabbit.PrintpartialF();
	rabbit.InitLocation();
	rabbit.InitSpeed();
	int count1 = 0;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	t2 = GetTickCount();
	printf("Use Time:%f\n", (t2 - t1)*1.0 / 1000);
	t1 = GetTickCount();
	do
		{
		rabbit.CalcNextLocation();
		rabbit.CalcNextForce();
		rabbit.CalcNextSpeed();
		rabbit.InitAgain();
		rabbit.CalcForce();
		rabbit.CalcpartialF();
		renderer->RemoveActor(actor);
		renderer->AddActor(actor);
		renderWindow->Render();
		ostringstream s1;
		s1 << ".\\Test\\test" << count1<<".vtk";
		string y = s1.str();
		cout << y << endl;
		writer->SetFileName(y.c_str());
		writer->SetInputData(rabbit.getPolydata());
		writer->Write();
		count1++;
		} while (count1<1000);
		t2 = GetTickCount();
		printf("Use Time:%f\n", (t2 - t1)*1.0 / 1000);

renderWindowInteractor->Start();
	
	}