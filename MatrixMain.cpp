#include <iostream>
#include <vector>
#include <string>
#include <ostream>
#include <math.h>
#include <algorithm>
#include "matrix.h"
#include <fstream>

class WrongSizeException : public exception {};
class VectorIsZeroException : public exception {};
class DetIsZeroException : public exception {};
class MatrixException : public exception {};

using namespace std;

int main() {

	/*
	vector< vector<double> > D{
		vector<double>{2, 4, -4},
		vector<double>{3, 8, 21},
		vector<double>{ -1, -2, -9}
	};
	vector< vector<double> > R{
		vector<double>{-7, 13, -9},
		vector<double>{ 2, 16, 6},
		vector<double>{ 10, 6, 0},
		vector<double>{ 5, 19, -1}
	};

	vector< vector<double> > R123{
	vector<double>{-7, 13, -9, 5},
	vector<double>{ 2, 16, 6, -3},
	};

	Matrix qwe(2, 4, R123);
	cout << " ~~~~~~~" << qwe.Rang() << endl;

	vector< vector<double> > V1{
		vector<double>{ 10, 11, 14}
	};
	//Matrix AV(1, 3, V1);
	vector< vector<double> > V2{
		vector<double>{0, 1, 4}
	};
	//Matrix BV(1, 3, V2);
	//Matrix CV = AV.VectorProduct(BV);	// Векторное произведение векторов
	//cout << "VectorProduct of AV and BV: " << CV << endl; //30 -40 10
	vector <vector<double> > Q{
		vector<double> {-5},
		vector<double> {3},
		vector<double> {7}
	};
	vector <vector<double> > W{
		vector<double> {3},
		vector<double> {6},
		vector<double> {0}
	};
	Matrix V11(3, 1, Q);
	Matrix V22(3, 1, W);
	Matrix Sc1(1, 3, V1);
	Matrix Sc2(1, 3, V2);
	Matrix OP = V11.OuterProduct(V22);
	Matrix DD(3, 3, D);


	cout << "Trace of DD: " << DD.Trace() << endl;
	cout << "Det of DD: " << DD.Determinant() << endl; // Вывожу матрицу после Гауса в этой же функции
	cout << "Scalar product of V11 and V22 (oba stolbtci): " << V11.ScalarProduct(V22) << endl;
	cout << "Scalar product of Sc1 and Sc2 (oba stroki): " << Sc1.ScalarProduct(Sc2) << endl;
	cout << "Outer product of V11 and V22: " << endl << OP; //Вектор-столбец на вектор-столбец (транспонирую второй)
	cout << "Euclid norm of Sc1 (stroka): " << Sc1.EuclidNormVector() << endl;
	cout << "Euclid norm of V11 (stolbetc): " << V11.EuclidNormVector() << endl;
	cout << "Max norm of Sc1: " << Sc1.MaxNormVector() << endl;
	cout << "Max norm of V11: " << V11.MaxNormVector() << endl;
	cout << "Matrix norm of DD: " << DD.EuclidNormMatrix() << endl;
	
	Matrix Rank(4, 3, R);
	Rank.Rang();

	cout << "Rang of Rank: " << Rank.Rang() << endl;
	cout << "Rang of DD: " << DD.Rang() << endl;

	Matrix ans = DD.InvertibleMatrix();
	cout << "Invertible matrix for DD:" << endl << ans << endl;
	Matrix qaz = ans * DD;
	cout << "Check for mistakes inv martix" << endl << qaz << endl;

	Matrix TranV11 = V11.Transpose();
	cout << "Traspose for V11: " << endl << TranV11 << endl;
	Matrix TranSc1 = Sc1.Transpose();
	cout << "Traspose for Sc1: " << endl << TranSc1 << endl;
	Matrix TranDD = DD.Transpose();
	cout << "Traspose for DD: " << endl << TranDD << endl;

	double AngleV11V22 = V11.Angle(V22);
	cout << "Angle between V11 and V22 (2 stolbca): " << endl << AngleV11V22 << endl;
	double AngleSc1Sc2 = Sc1.Angle(Sc2);
	cout << "Angle between Sc1 and Sc2 (2 stroki): " << endl << AngleSc1Sc2 << endl;
	*/
	/*
	ofstream in;
	in.open("in_file.txt");

	vector< vector<double> > D{
		vector<double>{2, 4, -4},
		vector<double>{3, 8, 21},
		vector<double>{ -1, -2, -9}
	};
	Matrix in_file(3, 3, D);

	if (in.is_open())
	{
		in << in_file << std::endl;
	}
	cout << in_file;
	

	ifstream out("out_of_file.txt");

	Matrix out_of_file(2, 2);
	out >> out_of_file;

	cout << out_of_file; // для корректной работы нужна пустая строка в конце текстового файла
	*/


	// 32 x 12

	Matrix data(32, 12);

	ifstream out_data("data.txt");

	out_data >> data;
	



	

	return 0;
}
