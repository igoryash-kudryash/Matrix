#pragma once

using namespace std;

class Matrix {
protected:
	vector< vector<double> > matrix;
	int n, m;										// n - строки, m - столбцы
public:
	Matrix();
	Matrix(int n, int m);
	Matrix(int n, int m, vector<vector<double>> X);
	~Matrix();
	int GetNSize();
	int GetMSize();
	double GetValue(int i, int j);
	void SetValue(int i, int j, double value);
	Matrix operator+ (Matrix& rvalue);				// Сложение
	Matrix operator* (double rvalue);				// Умножение на конст справа
	Matrix operator* (Matrix& rvalue);				// Умножение матриц
	Matrix Hadamar (Matrix& rvalue);				// Произведение Адамара
	Matrix operator- (Matrix& rvalue);				// Вычитание матриц
	double Trace();									// След матрицы
	double Determinant();							// Определитель
	double ScalarProduct(Matrix& value);			// Скалярное произведение векторов
	Matrix OuterProduct(Matrix& rvalue);			// Внешнее произведение векторов
	Matrix VectorProduct(Matrix& rvalue);			// Векторное произведение векторов (как площадь треугольника) ! NO NEED !
	double EuclidNormVector();						// Евклидова норма вектора
	double MaxNormVector();							// Максимальная норма вектора
	double EuclidNormMatrix();						// Евклидова норма матрицы
	double Angle(Matrix& value);					// Угол между векторами
	int Rang();										// Ранг матрицы
	Matrix InvertibleMatrix();						// Обратная матрица
	Matrix Transpose();								// Транспонирование матрицы

	void write_file(ofstream& o) const;
	void read_file(ifstream& i);



	//Matrix operator<<(string name);					// Вывод в текстовый
	//Matrix operator>>(string name);					// Ввод в текстовый
	//Matrix operator<(string name);					// Вывод в бинарный
	//Matrix operator>(string name);					// Ввод в бинарный
};

ostream& operator<<(ostream& out, Matrix& obj);		// Оператор вывода из текстового
istream& operator>>(istream& in, Matrix& obj);		// Оператор ввода в текстовый

Matrix operator* (double lvalue, Matrix rvalue);	// Умножение константы на матрицу

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

