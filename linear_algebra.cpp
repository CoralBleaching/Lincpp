#include<algorithm>
#include<numeric>
#include<functional>
#include<sstream>
#include<ostream>
#include<cmath>

#include "linear_algebra.hpp"

namespace alg {

	// VECTOR

	// Constructors

	Vector::Vector(size_t k, double val) : data_{ std::vector<double>(k, val) } {}
	Vector::Vector(const std::vector<double> v) : data_{ v } {}
	Vector::Vector(const std::initializer_list<double> l) : data_{ l } {}
	Vector::Vector(Iterator first, Iterator last) : data_{ std::vector<double>(first, last) } {}
	Vector::Vector(IteratorColumn first, IteratorColumn last) : data_{ std::vector<double>(first, last) } {}
	Vector::Vector(const Matrix& M) : data_{ M.getInternalStdVector() } {}

	// Shape methods

	size_t         Vector::size() const { return data_.size(); }
	bool           Vector::isScalar() const { return size() == 1; }
	void           Vector::clear() { data_.clear(); }

	Matrix Vector::row() const { return Matrix{ data_, 1, size() }; }
	Matrix row(Vector v) { return v.row(); }
	Matrix Vector::col() const { return Matrix{ data_, size(), 1 }; }
	Matrix col(Vector v) { return v.col(); }
	Matrix Vector::t() const { return row(); }
	inline Matrix t(Vector v) { return v.t(); }

	// Iterators and access

	std::vector<double>& Vector::getInternalStdVector() { return data_; }

	inline void Vector::push_back(double val) { data_.push_back(val); }

	void Vector::insert(Iterator it, double val)
	{
		std::ptrdiff_t dist = it - begin();
		std::vector<double>::iterator pos = data_.begin() + dist;
		data_.insert(pos, val);
	}

	void Vector::insert(Iterator it, const Vector& v)
	{
		std::ptrdiff_t dist = it - begin();
		std::vector<double>::iterator pos = data_.begin() + dist;
		std::copy(v.begin(), v.end(), std::inserter(data_, pos));
	}

	Vector Vector::concatenate(const Vector& v) const
	{
		Vector res(*this);
		res.insert(res.end(), v);
		return res;
	}

	Vector Vector::slice(size_t start, size_t finish)
	{
		size_t n = finish - start;
		if (n > size()) throw LinearAlgebraException("Slice exceeds Vector length.");
		Vector v(n);
		std::copy_n(begin() + start, n, v.begin());
		return v;
	}

	Iterator Vector::begin() const { return Iterator{ data_.begin()._Ptr }; }
	Iterator Vector::end() const { return Iterator{ data_.end()._Ptr }; }

	double& Vector::operator[](size_t i) { return data_[i]; }
	double& Vector::at(size_t i) { return data_.at(i); }

	inline std::string Vector::to_string() const { return toString(data_); }

	// Algebraic methods

	double Vector::norm(double power) const
	{
		return sqrt(std::accumulate<Iterator, double>(begin(), end(), 0, [=](double& accum, double& next)
			{
				return accum + pow(next, power);
			}));
	}

	// Operators

	// Unary operator

	Vector Vector::operator-() const
	{
		Vector v(size());
		std::transform(begin(), end(), v.begin(), std::negate<double>());
		return v;
	}

	// Operations with scalar

	Vector Vector::operator+(double t) const
	{
		Vector v{ *this };
		for (double& e : v) e += t;
		return v;
	}

	void Vector::operator+=(double t) { for (double& e : data_) e += t; }

	Vector operator+(double t, Vector v) { return v + t; }

	Vector Vector::operator-(double t) const
	{
		Vector v{ *this };
		for (double& e : v) e -= t;
		return v;
	}

	void Vector::operator-=(double t) { for (double& e : data_) e -= t; }

	Vector operator-(double t, Vector v)
	{
		for (double& e : v) e = t - e;
		return v;
	}

	Vector Vector::operator*(double t) const
	{
		Vector v{ *this };
		for (double& e : v) e *= t;
		return v;
	}

	void Vector::operator*=(double t) { for (double& e : data_) e *= t; }

	Vector operator*(double t, Vector v) { return v * t; }

	Vector Vector::operator/(double t) const
	{
		Vector v{ *this };
		for (double& e : v) e /= t;
		return v;
	}
	void Vector::operator/=(double t) { for (double& e : data_) e /= t; }

	Vector operator/(double t, Vector v) { for (double& e : v) e = t / e; return v; }

	// Operations with Vector

	double Vector::operator*(Vector v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}

	Vector Vector::operator+(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<double>());
		return v;
	}

	void Vector::operator+=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}

	Vector Vector::operator-(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<double>());
		return v;
	}

	void Vector::operator-=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}


	// Operations with Matrix


	Vector Vector::operator+(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<double>());
		return v;
	}


	void Vector::operator+=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}

	Vector operator+(const Matrix& M, Vector& v) { return v + M; }

	Vector Vector::operator-(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<double>());
		return v;
	}


	void Vector::operator-=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}

	Vector operator-(const Matrix& M, Vector v)
	{
		if (!M.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (v.size() != M.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		std::transform(v.begin(), v.end(), M.begin(), v.begin(), [](double a, double b) { return b - a; });
		return v;
	}

	Vector Vector::operator*(Matrix M) const
	{
		if (size() != M.nrows()) throw LinearAlgebraException("Can't perform dot product: Vector and Matrix have different number of columns.");
		Vector v;
		std::for_each(M.beginCol(), M.endCol(), [&](Column& col)
			{
				v.push_back(std::inner_product(col.begin(), col.end(), begin(), 0.));
			});
		return v;
	}

	Vector operator*(const Matrix& M, const Vector& v)
	{
		if (v.size() != M.ncols()) throw LinearAlgebraException("Can't perform dot product: Matrix and Vector have different number of lines.");
		Vector u;
		std::for_each(M.beginRow(), M.endRow(), [&](Row& row)
			{
				u.push_back(std::inner_product(row.begin(), row.end(), v.begin(), 0.));
			});
		return u;
	}

	// Other operations

	std::ostream& operator<<(std::ostream& ostream, Vector vector)
	{
		ostream << toString(vector);
		return ostream;
	}

	double norm(Vector v, double val) { return v.norm(val); }


	// MATRIX


	// Utils

	Matrix I(size_t m)
	{
		Matrix e{ m, m };
		for (size_t i = 0; i < m; i++) e[i][i] = (double)1;
		return e;
	}

	Vector arange(double end, double start, double step)
	{
		if ((start > end && step > 0.) || (start < end && step < 0.) || step == 0.)
			throw LinearAlgebraException("Invalid arguments to function range().");
		Vector rangelist(static_cast<size_t>(ceil((end - start)) / step), start);
		double loopstep = step;
		std::transform(rangelist.begin() + 1, rangelist.end(), rangelist.begin() + 1, [&](double& e) {
			double tmp = e + loopstep; loopstep += step; return tmp;
			});
		return rangelist;
	}

	inline typename Shape Matrix::shapeOf(const std::initializer_list<std::initializer_list<double>>& list)
	{
		return Shape{ list.size(), (*list.begin()).size() };
	}


	// Constructors

	alg::Matrix::Matrix(size_t k, alg::Shape shape) : data_{ std::vector<double>(k, 0) }, shape_{ shape } {}
	alg::Matrix::Matrix(size_t m, size_t n) : data_{ std::vector<double>(m * n, 0) }, shape_{ m, n } {}
	alg::Matrix::Matrix(Shape shape) : data_{ std::vector<double>(shape.m * shape.n, 0) }, shape_{ shape } {}
	alg::Matrix::Matrix(std::vector<double> data, Shape shape) : data_{ data }, shape_{ shape } {}
	alg::Matrix::Matrix(std::vector<double> data, size_t m, size_t n) : data_{ data }, shape_{ m, n } {
		if (m == -1) shape_ = { data.size(), 1 };
	}
	alg::Matrix::Matrix(alg::Vector vector, Shape shape) : data_{ vector.getInternalStdVector() }, shape_{ shape } {}
	alg::Matrix::Matrix(alg::Vector vector, size_t m, size_t n) : data_{ vector.getInternalStdVector() }, shape_{ m, n }
	{
		if (m == -1) shape_ = { vector.size(), 1 };
	}
	alg::Matrix::Matrix(const std::initializer_list<std::initializer_list<double>>& list)
	{
		setShape(shapeOf(list));
		for (auto& row : list)
			for (auto& elem : row)
				data_.emplace_back(elem);
	}

	// Shape methods


	inline size_t Matrix::size() const { return data_.size(); }
	inline Shape Matrix::getShape() const { return shape_; }
	inline size_t Matrix::nrows() const { return shape_.m; }
	inline size_t Matrix::ncols() const { return shape_.n; }
	inline void Matrix::setShape(size_t m, size_t n) { shape_.m = m; shape_.n = n; }
	inline void Matrix::setShape(Shape shape) { shape_ = shape; }
	Matrix Matrix::reshape(size_t m, size_t n) {
		//if (size() % m || size() % n) throw LinearAlgebraException("Cannot reshape matrix into desired shape.");
		return Matrix{ data_, { m, n} };
	};
	Matrix Matrix::reshape(Shape shape)
	{
		try { return reshape(shape.m, shape.n); }
		catch (LinearAlgebraException e) { throw e; }
	}
	inline bool Matrix::isScalar() const { return nrows() == 1 && ncols() == 1; }
	inline bool Matrix::isRow() const { return nrows() == 1; }
	inline bool Matrix::isColumn() const { return ncols() == 1; }
	inline bool Matrix::isVector() const { return isRow() || isColumn(); }
	inline void Matrix::clear() { data_.clear(); shape_ = { 0, 0 }; }

	void Matrix::push_back(const Vector& v)
	{
		if (v.size() != ncols()) throw LinearAlgebraException("Can't push back Vector of incompatible size.");
		setShape(nrows() + 1, ncols());
		for (double& e : v) data_.push_back(e);
	}

	Matrix Matrix::slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end) const
	{
		size_t m = row_end - row_begin;
		size_t n = column_end - column_begin;
		if (row_end > nrows() || column_end > ncols()) throw LinearAlgebraException("Slice exceeded matrix dimensions.");
		Matrix M = Matrix(m, n);
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				M[i][j] = data_[(i + row_begin) * ncols() + j + column_begin];
		return M;
	}


	void Matrix::insert(Iterator it, double val)
	{
		std::ptrdiff_t dist = it - begin();
		std::vector<double>::iterator idx = data_.begin() + dist;
		data_.insert(idx, val);
	}

	/*

	void Matrix::insertRows(IteratorRowVector& itIn_beg, IteratorRowVector& itOut_beg,
										IteratorRowVector& itOut_end)
	{
		if (itOut_beg.ncols() != ncols())
			throw LinearAlgebraException("Can't insert rows into matrix of different length.");
		auto it = itIn_beg.getIt();
		cout << "hello\n";
		auto dist = it - begin();
		cout << "hello " << dist << "\n";
		auto idx = matrix_.begin() + dist;
		cout << "hello " << *idx << "\n";
		std::for_each(itOut_beg, itOut_end, [&](auto& row)
			{
				cout << row << "\n";
				//cout << toString(matrix_) << "\n";
				cout << *row.begin() << " " << *(row.end() - 1) << "\n";
				//cout << toString(matrix_) << "\n";
				//matrix_.insert(idx, 69);
				//cout << toString(matrix_) << "\n";
				//matrix_.insert(idx, matrix_.begin(), matrix_.begin() + 2);
				//cout << toString(matrix_) << "\n";
				std::vector v;
				std::copy(row.begin(), row.end(), std::back_inserter(v));
				cout << toString(v) << "\n";
				matrix_.insert(idx, std::begin(v), std::end(v));
				cout << toString(matrix_) << "\n";
				//matrix_.insert(idx, row.begin(), row.end());
				idx += ncols();
			});
	}
	/**/
	/**

	void Matrix::insertColumns(IteratorColumnVector& itIn_beg, IteratorColumnVector& itOut_beg,
									IteratorColumnVector& itOut_end)
	{
		auto idx = itIn_beg.getIndex();
		auto it = matrix_.begin() + idx;
		std::for_each(itOut_beg, itOut_end, [&](auto& col)
			{
				for (size_t k = 0; k < col.size(); k++)
				{
					matrix_.insert(it + k * ncols(), col[k]);
				}
				it += 1;
			});
	}
	/**/

	Matrix Matrix::concatenate(const Matrix& other) const
	{
		Matrix mat{ size() + other.size(), { nrows(), ncols() + other.ncols() } };
		for (size_t i = 0; i < mat.nrows(); i++)
		{
			size_t j;
			for (j = 0; j < ncols(); j++)
			{
				mat[i][j] = data_[i * ncols() + j];
			}
			for (size_t k = 0; k < other.ncols(); k++, j++)
			{
				mat[i][j] = other.data_[i * other.ncols() + k];
			}
		}
		return mat;
	}

	// Iterators and access


	const double& Matrix::at(size_t i) const
	{
		if (i > size())
			throw LinearAlgebraException("Index exceeds flattened matrix length");
		return data_.at(i);
	}


	Iterator Matrix::begin() const { return Iterator{ data_.begin()._Ptr }; }


	Iterator Matrix::end() const { return Iterator{ data_.end()._Ptr }; }


	Row Matrix::row(size_t i) const
	{
		return Row{ begin() + i * ncols(), begin() + (i + 1) * ncols(), getShape() };
	}


	Column Matrix::col(size_t j) const
	{
		return Column
		{
			begin() + j,
			begin() + ncols() * nrows() + j,
			getShape()
		};
	}

	inline std::vector<double> Matrix::getInternalStdVector() const { return data_; }

	std::string Matrix::to_string() const
	{
		std::ostringstream ostr;
		ostr << "{";
		for (int i = 0; i < nrows(); i++)
		{
			ostr << ((i == 0) ? "{ " : " { ");
			for (int j = 0; j < ncols(); j++)
				ostr << at(j + i * ncols()) << " ";
			ostr << ((i == nrows() - 1) ? "}" : "}\n");
		}
		ostr << "}\n";
		return ostr.str();
	}


	IteratorRowVector Matrix::beginRow() const { return IteratorRowVector{ begin(), getShape() }; }


	IteratorRowVector Matrix::endRow() const { return IteratorRowVector{ begin() + nrows() * ncols(), getShape() }; }


	IteratorColumnVector Matrix::beginCol() const { return IteratorColumnVector{ this, 0, getShape() }; }


	IteratorColumnVector Matrix::endCol() const { return IteratorColumnVector{ &*this, ncols(), getShape() }; }

	// Algebraic methods
	/**/

	Matrix Matrix::inv() const
	{
		size_t m = nrows();
		size_t this_n = ncols();
		if (m != this_n) throw LinearAlgebraException("Matrix is not square.");
		Matrix M = Matrix(*this);
		M = M.concatenate(I(m));
		size_t n = M.ncols();
		for (size_t i = 0; i < m; i++)
		{
			size_t pivot = i;
			for (size_t k = i + 1; k < m; k++)
				if (fabs(M[k][i]) > fabs(M[i][i]))
					pivot = k;
			if (pivot > i)
			{
				for (size_t j = 0; j < n; j++)
				{
					double aux = M[i][j];
					M[i][j] = M[pivot][j];
					M[pivot][j] = aux;
				}
			}
			for (size_t k = i + 1; k < m; k++)
			{
				double mki = M[k][i] / M[i][i];
				//M[k][i] = 0;
				for (size_t j = i; j < n; j++)
				{
					M[k][j] -= mki * M[i][j];
				}
			}
		}
		for (int j = m - 1; j >= 0; j--)
		{
			double mjj = 1 / M[j][j];
			for (size_t k = j; k < n; k++)
				M[j][k] *= mjj;
			for (int i = j - 1; i >= 0; i--)
			{
				double mij = M[i][j];
				for (size_t k = j; k < n; k++)
				{
					double mij_mjk = -mij * M[j][k];
					M[i][k] -= mij * M[j][k];
				}
			}
		}
		M = M.slice(0, m, this_n, this_n * 2);
		return M;
	}


	double Matrix::determinant_recursion(Matrix M)
	{
		size_t m = M.nrows(), n = M.ncols();
		if (m != n)
			throw LinearAlgebraException("Matrix is not square.");
		if (n == 1)
			return M[0][0];
		if (n == 2)
			return  M[0][0] * M[1][1] - M[0][1] * M[1][0];
		else
		{
			double result = 0;
			for (size_t i = 0; i < n; i++)
			{
				Matrix left_submatrix = M.slice(1, n, 0, i);
				Matrix right_submatrix = M.slice(1, n, i + 1, n);
				Matrix submatrix = left_submatrix.concatenate(right_submatrix);
				result += std::pow(-1, i) * M[0][i] * determinant_recursion(submatrix);
			}
			return result;
		}
	}


	double Matrix::det() const { return determinant_recursion(*this); }


	Matrix Matrix::t() const
	{
		if (isVector())
		{
			if (isRow())
				return Matrix(data_, { ncols(), 1 }); // this constructor creates a column matrix.
			else
				if (isColumn())
				{
					return Matrix(data_, { 1, nrows() });
				}
		}
		else
		{
			size_t m = nrows();
			size_t n = ncols();
			Matrix trans{ n, m };
			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < m; j++)
					trans[i][j] = data_[j * ncols() + i];
			return trans;
		}
		throw LinearAlgebraException("doubleranspose error: matrix has zero dimensions.");
	}

	double Matrix::norm(double power = 2) const
	{
		return sqrt(std::accumulate<Iterator, double>(begin(), end(), 0, [&](double& accum, double& next)
			{
				return accum + pow(next, power);
			}));
	}

	// Operators

	// Unary operator


	Matrix Matrix::operator-() const
	{
		Matrix mat{ getShape() };
		std::transform(begin(), end(), mat.begin(), std::negate<double>());
		return mat;
	}

	// Operations with scalars


	Matrix Matrix::operator+ (double t) const
	{
		Matrix sum{ *this };
		for (double& e : sum) e += t;
		return sum;
	}


	void Matrix::operator+= (double t) { for (double& e : data_) e += t; }


	Matrix Matrix::operator- (double t) const
	{
		Matrix sum{ *this };
		for (double& e : sum) e -= t;
		return sum;
	}


	void Matrix::operator-= (double t) { for (double& e : data_) e -= t; }


	Matrix Matrix::operator* (double t) const
	{
		Matrix mul{ *this };
		for (double& e : mul) e *= t;
		return mul;
	}


	void Matrix::operator*= (double t) { for (double& e : data_) e *= t; }


	Matrix Matrix::operator/ (double t) const
	{
		Matrix mul{ *this };
		for (double& e : mul) e /= t;
		return mul;
	}


	void Matrix::operator/= (double t) { for (double& e : data_) e /= t; }

	// Operations with Matrix


	Matrix Matrix::operator+(const Matrix& other) const
	{
		if (nrows() != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
		if (ncols() != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
		Matrix sum{ getShape() };
		std::transform(begin(), end(), other.begin(), sum.begin(), std::plus<double>());
		return sum;
	}


	void Matrix::operator+=(const Matrix& other)
	{
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}


	Matrix Matrix::operator-(const Matrix& other) const
	{
		if (nrows() != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
		if (ncols() != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
		Matrix sum{ getShape() };
		std::transform(begin(), end(), other.begin(), sum.begin(), std::minus<double>());
		return sum;
	}


	void Matrix::operator-=(const Matrix& other)
	{
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}


	Matrix Matrix::operator*(const Matrix& other) const
	{
		size_t m = nrows();
		size_t n = ncols();
		size_t q = other.ncols();
		if (other.isScalar())
			return (*this) * other.data_[0];
		else if (isScalar())
			return other * data_[0];
		else if (n != other.nrows()) throw LinearAlgebraException("Multiplication error: matrices do not have appropriate dimensions.");
		else
		{
			//Matrix mul(m, q);
			Matrix mul(0, { m, q });
			std::for_each(beginRow(), endRow(), [&](Row& row)
				{
					std::for_each(other.beginCol(), other.endCol(), [&](Column& col)
						{
							mul.data_.push_back(
								std::inner_product(row.begin(), row.end(), col.begin(), 0.)
							);
						});
				});
			//for (size_t i = 0; i < m; i++)
			//	for (size_t j = 0; j < q; j++)
			//		for (size_t k = 0; k < n; k++)
			//			mul[i][j] += matrix_[i * ncols() + k] * other.matrix_[k * other.ncols() + j];
			return mul;
		}
	}


	void Matrix::operator*=(const Matrix& other)
	{
		size_t m = nrows();
		if (m != ncols())
			throw LinearAlgebraException("Multiplication error: assigning result of multiplication between non-square matrices to self matrix.");
		if (m != other.nrows() || m != other.ncols())
			throw LinearAlgebraException("Multiplication error: assining result of non-square matrix multiplication to left-hand-side matrix.");
		(*this) = (*this) * other;
	}

	// Other operations

	std::ostream& operator<<(std::ostream& ostream, Matrix matrix)
	{
		ostream << matrix.to_string(); return ostream;
	}

	/**/

	// ITERATORS

	// Row

		// Algebraic methods

	double Row::norm(double power) const
	{
		return sqrt(std::accumulate<Iterator, double>(begin(), end(), 0, [=](double& accum, double& next)
			{
				return accum + pow(next, power);
			}));
	}
	double norm(Row v, double val) { return v.norm(val); }

	// Operators

	// Unary operator


	Vector Row::operator-() const
	{
		Vector v(size());
		std::transform(begin(), end(), v.begin(), std::negate<double>());
		return v;
	}

	// Operations with scalar


	Vector Row::operator+(double t) const
	{
		Vector v(begin(), end());
		for (double& e : v) e += t;
		return v;
	}
	void Row::operator+=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e += t; });
	}
	Vector operator+(double t, Row v) { return v + t; }


	Vector Row::operator-(double t) const
	{
		Vector v(begin(), end());
		for (double& e : v) e -= t;
		return v;
	}
	void Row::operator-=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e -= t; });
	}
	Vector operator-(double t, Row v) { return v - t; }

	Vector Row::operator*(double t) const
	{
		Vector v(begin(), end());
		for (double& e : v) e *= t;
		return v;
	}

	void Row::operator*=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e *= t; });
	}
	Vector operator*(double t, Row v) { return v * t; }

	Vector Row::operator/(double t) const
	{
		Vector v(begin(), end());
		for (double& e : v) e /= t;
		return v;
	}
	void Row::operator/=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e /= t; });
	}
	Vector operator/(double t, Row v) { return v / t; }

	// Operations with Vector


	double Row::operator*(Vector v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}

	Vector Row::operator+(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<double>());
		return v;
	}


	void Row::operator+=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}


	Vector Row::operator-(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<double>());
		return v;
	}


	void Row::operator-=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}

	Vector operator-(Vector v, Row& r)
	{
		if (r.size() != v.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(v.begin(), v.end(), r.begin(), v.begin(), std::minus<double>());
		return v;
	}

	// Operations with Row


	double Row::operator*(Row v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}


	Vector Row::operator+(const Row& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<double>());
		return v;
	}


	void Row::operator+=(const Row& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}


	Vector Row::operator-(const Row& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<double>());
		return v;
	}


	void Row::operator-=(const Row& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}

	// Column

	// Algebraic methods

	double Column::norm(double power) const
	{
		return sqrt(std::accumulate<IteratorColumn, double>(begin(), end(), 0., [=](double& accum, double& next)
			{
				return accum + pow(next, power);
			}));
	}
	double norm(Column v, double val) { return v.norm(val); }

	// Operators

	// Unary operator


	Vector Column::operator-() const
	{
		Vector v(size());
		std::transform(begin(), end(), v.begin(), std::negate<double>());
		return v;
	}

	// Operations with scalar


	Vector Column::operator+(double t) const
	{
		Vector v{ begin(), end() };
		for (double& e : v) e += t;
		return v;
	}
	void Column::operator+=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e += t; });
	}
	Vector operator+(double t, Column v) { return v + t; }


	Vector Column::operator-(double t) const
	{
		Vector v{ begin(), end() };
		for (double& e : v) e -= t;
		return v;
	}
	void Column::operator-=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e -= t; });
	}
	Vector operator-(double t, Column v) { return v - t; }

	Vector Column::operator*(double t) const
	{
		Vector v{ begin(), end() };
		for (double& e : v) e *= t;
		return v;
	}

	void Column::operator*=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e *= t; });
	}
	Vector operator*(double t, Column v) { return v * t; }


	Vector Column::operator/(double t) const
	{
		Vector v{ begin(), end() };
		for (double& e : v) e /= t;
		return v;
	}
	void Column::operator/=(double t)
	{
		std::for_each(begin(), end(), [&](double& e) { e /= t; });
	}
	Vector operator/(double t, Column v) { return v / t; }

	// Operations with Vector


	double Column::operator*(Vector v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}


	Vector Column::operator+(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<double>());
		return v;
	}


	void Column::operator+=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}


	Vector Column::operator-(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<double>());
		return v;
	}


	void Column::operator-=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}


	// Operations with Row


	double Column::operator*(Column v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}


	Vector Column::operator+(const Column& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<double>());
		return v;
	}


	void Column::operator+=(const Column& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<double>());
	}


	Vector Column::operator-(const Column& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<double>());
		return v;
	}


	void Column::operator-=(const Column& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<double>());
	}

} // namespace::alg