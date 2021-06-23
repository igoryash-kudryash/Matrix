#include <iostream>
#include <vector>
#include <string>
#include <ostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "matrix.h"

class WrongSizeException : public exception {};
class VectorIsZeroException : public exception {};
class DetIsZeroException : public exception {};
class MatrixException : public exception {};

using namespace std;

double Round(double X) {
	double answer = X;
	if (abs(X - round(X)) < 0.000005) answer = round(X);
	return answer;
}

Matrix::Matrix(int n, int m) {
	this->n = n;
	this->m = m;
	for (int i = 0; i < n; i++) {
		matrix.push_back(vector<double>(m, 0));
	}
};
Matrix::Matrix() :Matrix(5, 5) {};
Matrix::Matrix(int n, int m, vector<vector<double> > X) {
	this->n = n;
	this->m = m;
	int checker1 = X[0].size();
	int checker2 = X.size();
	if (checker2 != n || checker1 != m) throw WrongSizeException();
		for (int i = 0; i < checker2; i++) {
			if (X[i].size() != m) throw WrongSizeException();
		}
	for (int i = 0; i < n; i++) {
		matrix.push_back(vector<double>(m, 0));
		for (int j = 0; j < m; j++) {
			matrix[i][j] = X[i][j];
		}
	}

};
Matrix::~Matrix() {
	matrix.clear();
};
IdentityMatrix::IdentityMatrix(int n): Matrix(n, n){
	for (int i = 0; i < n; i++) {
		matrix[i][i] = 1;
	}
};
DiagonalMatrix::DiagonalMatrix(int n, vector<double> X) : Matrix(n, n) {
	for (int i = 0; i < n; i++) {
		matrix[i][i] = X[i];
	}
};
UpperTrianMatrix::UpperTrianMatrix(int n, vector<double> X) : Matrix(n, n) {
	int k = 0;
	int size = 0;
	for (int i = 1; i <= n; i++) {
		size += i;
	}
	if (size != X.size()) throw WrongSizeException();
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix[i][j] = X[k];
			k++;
		}
	}
};
LowerTrianMatrix::LowerTrianMatrix(int n, vector<double> X) : Matrix(n, n) {
	int k = 0;
	int size = 0;
	for (int i = 1; i <= n; i++) {
		size += i;
	}
	if (size != X.size()) throw WrongSizeException();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= i; j++) {
			matrix[i][j] = X[k];
			k++;
		}
	}
};
SymmetricMatrix::SymmetricMatrix(int n, vector<double> X) : Matrix(n, n) {
	int size = 0;
	int k = 0;
	for (int i = 1; i <= n; i++) {
		size += i;
	}
	if (size != X.size()) throw WrongSizeException();
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix[i][j] = X[k];
			matrix[j][i] = X[k];
			k++;
		}
	}
};
double Matrix :: GetValue(int i, int j) {
	return matrix[i][j];
};
void Matrix::SetValue(int i, int j, double value) {
	matrix[i][j] = value;
}
int Matrix::GetNSize() {
	return n;
};
int Matrix::GetMSize() {
	return m;
};
Matrix Matrix::operator+(Matrix& rvalue) {
	Matrix answer(n, m);
	int checker = rvalue.GetNSize();
	int checker2 = rvalue.GetMSize();
	if (checker2 != m || checker != n) throw WrongSizeException();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			answer.SetValue(i, j, matrix[i][j] + rvalue.GetValue(i,j));
		}
	}
	return answer;
};
Matrix Matrix::operator*(double rvalue) {
	Matrix answer(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			answer.SetValue(i, j, matrix[i][j] * rvalue);
		}
	}
	return answer;
};
Matrix Matrix::operator*(Matrix& rvalue) {
	if (m != rvalue.GetNSize()) throw WrongSizeException();
	Matrix answer(n, rvalue.GetMSize());    // n*m X N*M
	double helper = 0;
	int k = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < rvalue.GetMSize(); j++) {
			for (int k = 0; k < m; k++) {
				helper = matrix[i][k] * rvalue.GetValue(k,j);
				answer.SetValue(i, j, answer.GetValue(i, j) + helper);
				helper = 0;
			}
		}
	}
	return answer;
};
Matrix Matrix::Hadamar(Matrix& rvalue) {
	Matrix answer(n, m);
	int checker = rvalue.GetNSize();
	int checker2 = rvalue.GetMSize();
	if (checker2 != m || checker != n) throw WrongSizeException();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			answer.SetValue(i, j, matrix[i][j] * rvalue.GetValue(i,j));
		}
	}
	return answer;
};
Matrix Matrix::operator-(Matrix& rvalue) {
	Matrix answer(n, m);
	int checker = rvalue.GetNSize();
	int checker2 = rvalue.GetMSize();
	if (checker2 != m || checker != n) throw WrongSizeException();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			answer.SetValue(i, j, matrix[i][j] - rvalue.GetValue(i, j));
		}
	}
	return answer;
};
Matrix operator* (double lvalue, Matrix rvalue) {
	int n = rvalue.GetNSize();
	int m = rvalue.GetMSize();
	Matrix answer(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			answer.SetValue(i, j, lvalue * rvalue.GetValue(i, j));
		}
	}
	return answer;
};
/*ostream& operator<<(ostream& out, Matrix& obj) {
	for (int i = 0; i < obj.GetNSize(); i++) {
		for (int j = 0; j < obj.GetMSize(); j++) {
			out << obj.GetValue(i, j) << '\t';
		}
		out << endl;
	}
	out << endl;
	return out;
};*/
double Matrix::Trace() {
	double answer = 0;
	for (int i = 0; i < n; i++) {
		answer += matrix[i][i];
	}
	return answer;
};
int Gauss(Matrix& a) {	// реализовано как отдельная функция, а не как метод, чтобы удобно было вызывать в других методах
						// меняю переданную матрицу, а возвращаю количеству раз, когда мы меняли строки местами (для определителя)
	int n = a.GetNSize();						
	int m = a.GetMSize();							
	int k = 0, index;
	int counter = 0;
	double max, helper, h1;
	while (k < n) {
		if (k < n && k < m) max = abs(a.GetValue(k, k));
		else break;
		index = k;
		for (int i = k + 1; i < n; i++) { // поиск строки с максимальным первым элементом
			if (abs(a.GetValue(i, k)) > max) {
				max = a.GetValue(i, k);
				index = i;
			}
		}

		for (int j = 0; j < m; j++) { // перестановка строки с максимальным наверх
			double tmp1 = a.GetValue(k, j);
			a.SetValue(k, j, a.GetValue(index, j));
			a.SetValue(index, j, tmp1);

			if (index != k) counter++;
		}

		for (int i = k + 1; i < n; i++) {
			if (a.GetValue(i, k) != 0) {
				helper = a.GetValue(i, k) / max;
				for (int j = k; j < m; j++) {
					h1 = a.GetValue(i, j) - helper * a.GetValue(k, j);
					a.SetValue(i, j, Round(h1));
				}
			}
		}
		k++;
	}
	return pow(-1, counter);
}
double Matrix::Determinant() {
	Matrix t(*this);
	if (n != m) throw WrongSizeException();
	int sign = Gauss(t);
	double answer = 1;
	//cout << "Our matrix after Gauss : " << endl << t;
	for (int i = 0; i < n; i++) {
		answer *= t.GetValue(i, i);
	}
	return answer*sign;
}
double Matrix::ScalarProduct(Matrix& value) { //неважно вектор-строки или вектор-столбцы
	double answer = 0;
	if (value.GetNSize() == 1 && n == 1) {
		if (value.GetMSize() == m) {
			for (int i = 0; i < m; i++) {
				answer += matrix[0][i] * value.GetValue(0, i);
			}
		}
		else { throw WrongSizeException(); }
	}
	else {
		if (value.GetMSize() == 1 && m == 1) {
			if (value.GetNSize() == n) {
				for (int i = 0; i < n; i++) {
					answer += matrix[i][0] * value.GetValue(i, 0);
				}
			}
			else { throw WrongSizeException(); }
		}
		else { throw WrongSizeException(); }
	}

	return answer;
}
Matrix Matrix::OuterProduct(Matrix& rvalue) { // AV = this 1x4  // rvalue=BV 1x4
	Matrix answer;
	Matrix helper = rvalue.Transpose(); //helper = 4x1 BV
	if (this->GetMSize() != 1 || helper.GetNSize() != 1 || (this->GetNSize() != helper.GetMSize())) throw WrongSizeException();
	answer = (*this) * helper;
	return answer;
}
Matrix Matrix::VectorProduct(Matrix& value) { // NO NEED
	if (this->GetNSize() != 1 || value.GetNSize() != 1 || (this->GetMSize() != value.GetMSize())) throw WrongSizeException();
	Matrix answer(1,m);
	for (int i = 0; i < m; i++) {
		answer.SetValue(0, i, this->GetValue(0, (i+1) % m) * value.GetValue(0, (i+2) % m) 
			- this->GetValue(0, (i+2) % m) * value.GetValue(0, (i+1) % m));
	}
	return answer;
} // NO NEED !
double Matrix::EuclidNormVector() { // неважно, вектор в столбец или в строку
	double answer = 0;
	if (n == 1) {
		for (int i = 0; i < m; i++) {
			answer += pow(matrix[0][i], 2);
		}
	}
	else {
		if (m == 1 ) {
			for (int i = 0; i < n; i++) {
				answer += pow(matrix[i][0], 2);
			}
		}
		else { throw WrongSizeException(); }
	}
	return pow(answer,0.5);
}
double Matrix::MaxNormVector() { 
	double max = abs(matrix[0][0]);
	if (n == 1) {
		for (int i = 1; i < m; i++) {
			if (abs(matrix[0][i]) > max) max = abs(matrix[0][i]);
		}
	}
	else {
		if (m == 1) {
			for (int i = 1; i < n; i++) {
				if (abs(matrix[i][0]) > max) max = abs(matrix[i][0]);
			}
		}
		else { throw WrongSizeException(); }
	}
	return max;
}
double Matrix::EuclidNormMatrix() {
	double answer = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			answer += pow(matrix[i][j], 2);
		}
	}
	return pow(answer,0.5);
}
double Matrix::Angle(Matrix& value) { // без разницы вектор-строки или вектор-столбцы
	double Scalar, l1, l2, answer;
	if (n == 1 && value.GetNSize() == 1) {
		if (m == value.GetMSize()) {
			l1 = this->EuclidNormVector();
			l2 = value.EuclidNormVector();
			if (l1 == 0 || l2 == 0) throw VectorIsZeroException();
			Scalar = this->ScalarProduct(value);
			answer = (Scalar / (l1 * l2));
			return answer;
		}
		else throw WrongSizeException();
	}
	if (m == 1 && value.GetMSize() == 1) {
		if (n == value.GetNSize()) {
			l1 = this->EuclidNormVector();
			l2 = value.EuclidNormVector();
			if (l1 == 0 || l2 == 0) throw VectorIsZeroException();
			Scalar = this->ScalarProduct(value);
			answer = acos(Scalar / (l1 * l2));
			return answer;
		}
		else throw WrongSizeException();
	}
	throw WrongSizeException();
}
int Matrix::Rang() {
	Matrix t(*this);
	Gauss(t);
	int n = t.GetNSize(), m = t.GetMSize();
	int counter = 0;
	int rank = n;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (t.GetValue(i, j) != 0) break;
			else {
				for (int k = 0; k < m; k++) {
					if (t.GetValue(i, k) == 0) continue;
					else counter = 1;
				}
				if (counter == 0) {
					rank--;
					break;
				}
				counter = 0;
			}
		}
	}
	return rank;
}
void Swap(Matrix& A, int x, int y) { // используется только для обратной матрицы
	double  t = 0;
	int c = A.GetMSize();
	for (int i = 0; i < c; ++i) {
		t = A.GetValue(x, i);
		A.SetValue(x, i, A.GetValue(y, i));
		A.SetValue(y, i, t);
	}
}
Matrix Matrix::InvertibleMatrix() {
	Matrix t(*this);
	if (t.Determinant() == 0) throw DetIsZeroException();
	int n = t.GetNSize(), col = 0;
	double k, max;
	IdentityMatrix inv(n);

	for (int i = 0; i + col < n; ++i) {
		if (!t.GetValue(i, i + col)) {
			for (int j = i + 1; j < n; ++j) {
				if (t.GetValue(j, i + col)) {
					Swap(t, i, j);
					Swap(inv, i, j);
					break;
				}
			}
			if (!t.GetValue(i, i + col)) {
				++col;
				--i;
			}
		}
	}

	for (int i = 0; i < n - 1; ++i) {
		max = t.GetValue(i, i);
		if (max) { 
			for (int j = i + 1; j < n; ++j) {
				k = -t.GetValue(j, i) / max;
				for (int g = 0; g < n; ++g) {
					t.SetValue(j, g, t.GetValue(j, g) + k * t.GetValue(i, g));
					inv.SetValue(j, g, inv.GetValue(j, g) + k * inv.GetValue(i, g));
				}
			}
		} 
	}

	for (int i = n - 1; i > 0; --i) {
		max = t.GetValue(i, i);
		for (int j = i - 1; j >= 0; --j) {
			k = -t.GetValue(j, i) / max;
			for (int g = 0; g < n; ++g) {
				t.SetValue(j, g, t.GetValue(j, g) + k * t.GetValue(i, g));
				inv.SetValue(j, g, inv.GetValue(j, g) + k * inv.GetValue(i, g));
			}
		}
	}

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			inv.SetValue(i, j, inv.GetValue(i, j) / t.GetValue(i, i));

	return inv;
}
Matrix Matrix::Transpose() {
	Matrix helper(m, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			helper.SetValue(j, i, matrix[i][j]);
		}
	}
	return helper;
};

/*
ostream& operator<<(ostream& out, Matrix& obj) {
	const int k = 10000;

	for (int i = 0; i < obj.GetNSize(); i++) {
		for (int j = 0; j < obj.GetMSize() - 1; j++)
			out << round(obj.GetValue(i, j) * k) / k << '\t';
		out << round(obj.GetValue(i, obj.GetMSize() - 1) * k) / k;
		out << endl;
	}
	out << endl;

	return out;
}
istream& operator>>(istream& i, Matrix& obj) {
	vector<double> helper;
	vector<vector<double> > ans;
	double x;
	int k;
	string line;
	stringstream ss;

	getline(i, line);
	int m = count(line.begin(), line.end(), '\t') + 1, n = 0;

	while (line != "") {
		++n;

		ss << line;
		for (k = 0; k < m; ++k) {
			ss >> x;
			helper.push_back(x);
		}

		getline(i, line);
		ans.push_back(helper);
		helper.clear();
		ss.str(string());
		ss.clear();
	}

	obj = Matrix(n, m, ans);

	return i;
}
void Matrix::WriteFile(ofstream& out) {
	out.write((char*)&n, sizeof(int));
	out.write((char*)&m, sizeof(int));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			out.write((char*)&(matrix[i][j]), sizeof(double));
		}
	}
}
void Matrix::ReadFile(ifstream& in) {
	double x;

	in.read((char*)&n, sizeof(int));
	in.read((char*)&m, sizeof(int));

	vector <double> helper;
	vector <vector<double> > ans;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			in.read((char*)&x, sizeof(double));
			helper.push_back(x);
		}
		ans.push_back(helper);
		helper.clear();
	}

	*this = Matrix(n, m, ans);
}
*/

/*
Matrix Matrix::operator<<(string name) {
	ifstream in(name);
	vector<vector<double> > ans;
	vector<double> X;
	double helper;

	if (in.is_open()) {
		int count = 0;
		int temp;
		int k = 1, i = 0, pos = 0;
		string s;

		int count_space = 0, count_space1 = 0;
		char symbol, symbol1;

		while (!in.eof()) { // Считаем строки и столбцы отдельно и получаем размер матрицы
			in.get(symbol);
			if (symbol == '\t') count_space++;
			if (symbol == '\n') k++;
		}
		in.close();
		ifstream in1(name);

		while (!in1.eof()) {
			in1 >> temp;
			count++;
		}
		in1.close();

		if (count != count_space + k) { throw MatrixException(); }
		ifstream in2(name);

		m = (count_space + k) / k;
		n = k;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				in2 >> helper;
				X.push_back(helper);
			}
			ans.push_back(X);
			X.clear();
		}
		return Matrix(n, m, ans);
		in2.close();
	}
	else {
		throw MatrixException();
		return Matrix(n, m);
	}
}
Matrix Matrix::operator>>(string name) {
	fstream inOut;
	vector<vector<double> > ans;
	vector<double> X;
	Matrix t(*this);

	inOut.open(name, ios::out);
	if (inOut.is_open()){
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				inOut << t.GetValue(i, j) << '\t';
			}
			inOut << endl;
		}
		inOut.close();		//под конец закроем файла
		return Matrix(n, m, matrix);
	}
	else {   //Если открытие файла прошло не успешно
		throw MatrixException();
		return Matrix(n, m, matrix);
	}
}
Matrix Matrix::operator<(string name) {
	ifstream inOut(name, ios::binary);
	vector<vector<double> > ans;
	vector<double> X;

	if (inOut.is_open()){
		double x;
		inOut.read((char*)&n, sizeof(int));
		inOut.read((char*)&m, sizeof(int));

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				inOut.read((char*)&x, sizeof(double));
				X.push_back(x);
			}
			ans.push_back(X);
			X.clear();
		}

		inOut.close();
		return Matrix(n, m, ans);
	}
	else{	//Если открытие файла прошло не успешно
		throw MatrixException();
		return Matrix(n, m, ans);
	}
}
Matrix Matrix::operator>(string name) {
	ofstream inOut(name, ios::out || ios::binary);
	double x;

	if (inOut.is_open()) {
		inOut.write((char*)&n, sizeof(int));
		inOut.write((char*)&m, sizeof(int));

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				x = matrix[i][j];
				inOut.write((char*)&(x), sizeof(double));
			}
		}
		inOut.close();
	}
	else {	//Если открытие файла прошло не успешно
		throw MatrixException();
		return Matrix(n, m, matrix);
	}
}
*/

/*Matrix operator*( Matrix& lvalue, Matrix& rvalue) {
	if (lvalue.GetMSize() != rvalue.GetNSize()) throw WrongSizeException();
	Matrix answer(lvalue.GetNSize(), rvalue.GetMSize());    // n*m X N*M
	double helper = 0;
	int k = 0;
	for (int i = 0; i < lvalue.GetNSize(); i++) {
		for (int j = 0; j < rvalue.GetMSize(); j++) {
			for (int k = 0; k < lvalue.GetMSize(); k++) {
				helper = lvalue.GetValue(i, k) * rvalue.GetValue(k, j);
				answer.SetValue(i, j, answer.GetValue(i, j) + helper);
				helper = 0;
			}
		}
	}
	return answer;
}*/




ostream& operator<<(ostream& o, Matrix& obj) {
	const int k = 10000;

	for (int i = 0; i < obj.GetNSize(); i++) {
		for (int j = 0; j < obj.GetMSize() - 1; j++)
			o << round(obj.GetValue(i, j) * k) / k << '\t';
		o << round(obj.GetValue(i, obj.GetMSize() - 1) * k) / k;
		o << endl;
	}
	o << endl;

	return o;
}

istream& operator>>(istream& i, Matrix& obj) {
	vector<double> t;
	double x;
	int k;
	string line;
	stringstream ss;

	getline(i, line);
	int c = count(line.begin(), line.end(), '\t') + 1, r = 0;
	
	int counter = 0;

	while (line != "") {

		
		++r;

		ss << line;
		for (k = 0; k < c; ++k) {
			ss >> x;
			t.push_back(x);
		}

		getline(i, line);
		ss.str(string());
		ss.clear();
	}

	vector<double> vect;
	vector<vector<double> > helper;

	k = 0;
	for (int j = 0; j < r; j++) {
		for (int i = 0; i < c; i++) {
			vect.push_back(t[k]);
			k++;
		}
		helper.push_back(vect);
		vect.clear();
	}

	obj = Matrix(r, c, helper);

	return i;
}

void Matrix::write_file(ofstream& o) const {
	o.write((char*)&n, sizeof(int));
	o.write((char*)&m, sizeof(int));

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			o.write((char*)&(matrix[i][j]), sizeof(double));
}

void Matrix::read_file(ifstream& i) {
	double x;

	i.read((char*)&n, sizeof(int));
	i.read((char*)&m, sizeof(int));

	vector<double> vect;
	vector<vector<double> >helper;

	for (int k = 0; k < n; k++) {
		for (int j = 0; j < m; j++) {
			i.read((char*)&x, sizeof(double));
			vect.push_back(x);
		}
		helper.push_back(vect);
		vect.clear();
	}

	*this = Matrix(n, m, helper);
}

PCA::PCA(string name) {
	ifstream in(name);
	in >> info;
	in.close();
}

PCA::PCA() {}

PCA::~PCA() {}

void PCA::autoscaling() {
	int r = info.GetNSize(), c = info.GetMSize();
	double m, s;

	for (int j = 0; j < c; ++j)
	{
		m = 0;
		for (int i = 0; i < r; ++i)
			m += info.GetValue(i, j);
		m /= r;

		s = 0;
		for (int i = 0; i < r; ++i)
			s += (info.GetValue(i, j) - m) * (info.GetValue(i, j) - m);
		s /= r - 1;
		s = sqrt(s);

		for (int i = 0; i < r; ++i)
			info.SetValue(i, j, (info.GetValue(i, j) - m) / s);
	}
}

void PCA::NIPALS(Matrix& T, Matrix& P, Matrix& E, int A) {
	autoscaling();
	//cout << info;

	int r = info.GetNSize(), c = info.GetMSize();
	vector<vector<double> > ans;
	vector<double>helper;
	const double threshold = pow(0.1, 5);
	Matrix E_a(info);
	int a = A ? A : c;

	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; j++){
			helper.push_back(E_a.GetValue(i, 0));
		}
		ans.push_back(helper);
		helper.clear();
	}

	vector<vector<double> > T_ans;
	vector<double> T_helper;
	vector<vector<double> > P_ans;
	vector<double> P_helper;

	Matrix t(r, 1, ans), p, t_old;

	int k = 0;
	while (k < a) {
		t_old = t;
		p = (t.Transpose() * E_a * (1 / t.ScalarProduct(t))).Transpose();
		p = p * sqrt(1 / p.ScalarProduct(p));
		t = E_a * p * (1 / p.ScalarProduct(p));

		//https://folk.uio.no/henninri/pca_module/pca_nipals.pdf
		if (abs(t_old.ScalarProduct(t_old) - t.ScalarProduct(t)) < threshold * t.ScalarProduct(t)) {

			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < a; j++) {
					T_helper.push_back(t.GetValue(i, 0));
				}
				T_ans.push_back(T_helper);
				T_helper.clear();
			}
			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < a; j++) {
					P_helper.push_back(t.GetValue(i, 0));
				}
				P_ans.push_back(P_helper);
				P_helper.clear();
			}
			++k;
			Matrix p_T_helper = p.Transpose();
			p_T_helper = t * p_T_helper;
			E_a = E_a - p_T_helper;
		}
	}

	T = Matrix(r, a, T_ans);
	P = Matrix(c, a, P_ans);
	Matrix P_T_helper = P.Transpose();
	P_T_helper = T * P_T_helper;
	E = info - P_T_helper;
}

Matrix PCA::leverage(Matrix& T) {
	int r = T.GetNSize(), c = T.GetMSize();
	vector<double>vect (c);	
	vector<double>H_vect (r);
	Matrix t;
	vector<vector<double> > helper;
	vector<vector<double> > helper2;
	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			vect[j] = T.GetValue(i, j);
		}

		helper.push_back(vect);

		t = Matrix(c, 1 , helper);

		Matrix h1 = (T.Transpose() * T);
		Matrix h2 = t.Transpose() * h1;
		H_vect[i] = (h2.InvertibleMatrix() * t).GetValue(0, 0);

		helper.clear();
	}
	helper2.push_back(H_vect);

	return Matrix(r, 1, helper2);
}

Matrix PCA::variation(Matrix& E) {
	int r = info.GetNSize(), c = info.GetMSize();
	vector<double> V_vect (r);
	vector< vector<double> > helper;

	for (int i = 0; i < r; ++i) {
		V_vect[i] = 0;
		for (int j = 0; j < c; ++j)
			V_vect[i] += E.GetValue(i, j) * E.GetValue(i, j);
	}

	helper.push_back(V_vect);
	return Matrix(r, 1, helper);
}

double PCA::TRV(Matrix& E) {
	int r = info.GetNSize(), c = info.GetMSize();
	Matrix v = variation(E);
	double s = 0;

	for (int i = 0; i < r; ++i)
		s += v.GetValue(i, 0);
	s /= r;
	s /= c;

	return s;
}

double PCA::ERV(Matrix& E)  {
	return 1 - TRV(E) / TRV(info);
}