#pragma once

using namespace std;

class Matrix {
protected:
	vector< vector<double> > matrix;
	int n, m;										// n - ������, m - �������
public:
	Matrix();
	Matrix(int n, int m);
	Matrix(int n, int m, vector<vector<double>> X);
	~Matrix();
	int GetNSize();
	int GetMSize();
	double GetValue(int i, int j);
	void SetValue(int i, int j, double value);
	Matrix operator+ (Matrix& rvalue);				// ��������
	Matrix operator* (double rvalue);				// ��������� �� ����� ������
	Matrix operator* (Matrix& rvalue);				// ��������� ������
	Matrix Hadamar (Matrix& rvalue);				// ������������ �������
	Matrix operator- (Matrix& rvalue);				// ��������� ������
	double Trace();									// ���� �������
	double Determinant();							// ������������
	double ScalarProduct(Matrix& value);			// ��������� ������������ ��������
	Matrix OuterProduct(Matrix& rvalue);			// ������� ������������ ��������
	Matrix VectorProduct(Matrix& rvalue);			// ��������� ������������ �������� (��� ������� ������������) ! NO NEED !
	double EuclidNormVector();						// ��������� ����� �������
	double MaxNormVector();							// ������������ ����� �������
	double EuclidNormMatrix();						// ��������� ����� �������
	double Angle(Matrix& value);					// ���� ����� ���������
	int Rang();										// ���� �������
	Matrix InvertibleMatrix();						// �������� �������
	Matrix Transpose();								// ���������������� �������

	void write_file(ofstream& o) const;
	void read_file(ifstream& i);



	//Matrix operator<<(string name);					// ����� � ���������
	//Matrix operator>>(string name);					// ���� � ���������
	//Matrix operator<(string name);					// ����� � ��������
	//Matrix operator>(string name);					// ���� � ��������
};

ostream& operator<<(ostream& out, Matrix& obj);		// �������� ������ �� ����������
istream& operator>>(istream& in, Matrix& obj);		// �������� ����� � ���������

Matrix operator* (double lvalue, Matrix rvalue);	// ��������� ��������� �� �������

//Matrix operator* (Matrix& lvalue, Matrix& rvalue);

class IdentityMatrix : public Matrix {
public:
	IdentityMatrix(int n);
};

class DiagonalMatrix : public Matrix {
public:
	DiagonalMatrix(int n, vector<double> X);
};

class UpperTrianMatrix : public Matrix {
public:
	UpperTrianMatrix(int n, vector<double> X);
};

class LowerTrianMatrix : public Matrix {
public:
	LowerTrianMatrix(int n, vector<double> X);
};

class SymmetricMatrix : public Matrix {
public:
	SymmetricMatrix(int n, vector<double> X);
};


class PCA {
public:
	PCA(std::string);
	PCA();
	~PCA();

	void autoscaling();
	void NIPALS(Matrix&, Matrix&, Matrix&, int = 0);
	Matrix leverage(Matrix&);
	Matrix variation(Matrix&);
	double TRV(Matrix&);
	double ERV(Matrix&);
private:
	Matrix info;
};

