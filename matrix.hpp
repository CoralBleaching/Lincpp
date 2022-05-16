#pragma once

#include<vector>
#include<algorithm>
#include<numeric>
#include<functional>
#include<sstream>
#include<ostream>
#include<cmath>

#include "iterators.hpp"

namespace alg {

	template<class T> class Vector;
	template<class T> Matrix<T> I(size_t);
	class LinearAlgebraException;
	template<typename Container> std::string toString(Container con);
	template<class T, class... Args> Vector<T> range(Args... args);

	template<class T>
	class Matrix
	{

		// DECLARATIONS
	private:

		static T determinant_recursion(Matrix<T>);

	protected:

		std::vector<T> matrix_;
		Shape shape_;

		static inline Shape shapeOf(const std::initializer_list<std::initializer_list<T>>&);

	public: // PRELOAD ALL ROWS AND COLUMNS FOR SPEED EFFICIENCY?
		Matrix(size_t k = 0, Shape shape = { 0, 0 }) : matrix_{ std::vector<T>(k, 0) }, shape_{ shape } {}
		Matrix(size_t m, size_t n) : matrix_{ std::vector<T>(m * n, 0) }, shape_{ m, n }  {}
		Matrix(Shape shape) : matrix_{ std::vector<T>(shape.m * shape.n, 0) }, shape_{ shape } {}
		Matrix(std::vector<T> matrix, Shape shape = { matrix.size(), 1 }) : matrix_{ matrix }, shape_{ shape } {}
		Matrix(std::vector<T> matrix, size_t m = matrix.size(), size_t n = 1) : matrix_{ matrix }, shape_{ m, n } {}
		Matrix(Vector<T> vector, Shape shape = { vector.size(), 1 }) : matrix_{ vector.getInternalStdVector() }, shape_{ shape } {}
		Matrix(Vector<T> vector, size_t m = vector.size(), size_t n = 1) : matrix_{ vector.getInternalStdVector() }, shape_{ m, n } {}
		Matrix(const std::initializer_list<std::initializer_list<T>>&);
		Matrix(const Matrix&);

		//friend std::ostream& operator<< (std::ostream&, Matrix<T>&);

		// Shape methods

		inline size_t size() const { return matrix_.size(); }
		inline Shape getShape() const { return shape_; }
		inline size_t nrows() const { return shape_.m; }
		inline size_t ncols() const { return shape_.n; }
		inline void setShape(size_t m, size_t n) { shape_.m = m; shape_.n = n; }
		inline void setShape(Shape shape) { shape_ = shape; }
		Matrix<T> reshape(size_t m, size_t n);
		Matrix<T> reshape(Shape shape);
		inline bool isScalar() const { return nrows() == 1 && ncols() == 1; }
		inline bool isRow() const { return nrows() == 1; }
		inline bool isColumn() const { return ncols() == 1; }
		inline bool isVector() const { return isRow() || isColumn(); }
		inline void clear() { matrix_.clear(); shape_ = { 0, 0 }; }

		//class ReferenceMatrix;

		void push_back(const Vector<T>& v)
		{
			if (v.size() != ncols()) throw LinearAlgebraException("Can't push back Vector of incompatible size.");
			setShape(nrows() + 1, ncols());
			for (T& e : v) matrix_.push_back(e);
		}
		void insert(Iterator<T> it, const T& val);
		//void insertRows(IteratorRowVector& itIn_beg, IteratorRowVector& itOut_beg, IteratorRowVector& itOut_end);
		//void insertColumns(IteratorColumnVector& itIn_beg, IteratorColumnVector& itOut_beg, IteratorColumnVector& itOut_end);
		Matrix<T> concatenate(const Matrix<T>&) const;
		//Matrix<T> concatenateBelow(const Matrix<T>&);
		Matrix<T> slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end) const;
		//virtual Matrix<T> slice(const arrayIterator<T>& it_b, const arrayIterator<T>& it_e) const;

		// Iterators and access

		inline Iterator<T> begin() const;
		inline Iterator<T> end() const;
		inline IteratorRowVector<T> beginRow() const;
		inline IteratorRowVector<T> endRow() const;
		inline IteratorColumnVector<T> beginCol() const;
		inline IteratorColumnVector<T> endCol() const;

		T& at(size_t i);
		//Row operator[](size_t i) const;
		T* operator[](size_t i) { return matrix_.data() + i * ncols(); }
		T& operator()(size_t i, size_t j) { return matrix_[i * ncols() + j]; }
		// ReferenceMatrix operator[](std::initializer_list<size_t>, std::initializer_list<size_t>) const;
		Row<T> row(size_t i) const;
		Column<T> col(size_t j) const;

		inline std::vector<T> getInternalStdVector() const { return matrix_; }
		/**/
	public:

		// Algebraic methods

		Matrix<T> inv() const;
		friend Matrix<T> inv(Matrix<T> M) { return M.inv(); }
		T det() const;
		friend T det(Matrix<T> M) { return M.det(); }
		Matrix<T> t() const;
		friend Matrix<T> t(Matrix<T> M) { return M.t(); }
		T norm(T) const;
		friend T norm(Matrix<T> M, T pow = 2) { return M.norm(pow); }

		// Operators

		// Unary operator

		Matrix<T> operator- () const;

		// Operations with scalars:

		Matrix<T> operator+ (T) const;
		Matrix<T> operator- (T) const;
		Matrix<T> operator* (T) const;
		Matrix<T> operator/ (T) const;
		void operator+= (T);
		void operator-= (T);
		void operator*= (T);
		void operator/= (T);

		friend Matrix<T> operator+ (T t, const Matrix<T>& M) { return M + t; }
		friend Matrix<T> operator- (T t, const Matrix<T>& M) { return M - t; }
		friend Matrix<T> operator* (T t, const Matrix<T>& M) { return M * t; }
		friend Matrix<T> operator/ (T t, const Matrix<T>& M) { return M / t; }

	public:

		// Operations with Matrix

		Matrix<T> operator+ (const Matrix<T>&) const;
		void operator+= (const Matrix<T>&);
		Matrix<T> operator- (const Matrix<T>&) const;
		void operator-= (const Matrix<T>&);
		Matrix<T> operator* (const Matrix<T>&) const;
		void operator*= (const Matrix<T>&);
	};

	/**/
	
	// DEFINITIONS

	// Utils

	class LinearAlgebraException : public std::exception
	{
	public:

		const char* message_;
		LinearAlgebraException(const char* message = "") { message_ = message; }
		virtual const char* what() const throw() { return message_; }
	};

	template<class T>
	Matrix<T> I(size_t m)
	{
		Matrix<T> e{ m, m };
		for (size_t i = 0; i < m; i++) e[i][i] = (T)1;
		return e;
	}

	template<typename Container>
	std::string toString(Container con)
	{
		std::ostringstream str;
		str << "{ ";
		for (auto it = con.begin(); it != con.end(); it++) str << *it << " ";
		str << "}";
		return str.str();
	}

	template<class T, class... Args>
	Vector<T> range(Args... args)
	{
		std::vector<T> arglist{ args... };
		T start, end, step = 1;
		int numargs = arglist.size();
		switch (numargs)
		{
		case 3:
			step = arglist[2];
		case 2:
			end = arglist[1]; start = arglist[0];
			break;
		case 1:
			end = arglist[0]; start = 0;
			break;
		default:
			throw LinearAlgebraException("range() needs to be called with at least 1 argument, and maximum 3.");
		}
		if ((start > end && step > 0) || (start < end && step < 0) || step == 0)
			throw LinearAlgebraException("Invalid arguments to function range().");
		Vector<T> rangelist(ceil((end - start) / (float)step), start);
		T loopstep = step;
		std::transform(rangelist.begin() + 1, rangelist.end(), rangelist.begin() + 1, [&](T& e) {
			T tmp = e + loopstep; loopstep += step; return tmp;
			});
		return rangelist;
	}

	template<class T>
	std::ostream& operator<< (std::ostream& output, Matrix<T>& M)
	{
		std::ostringstream str;
		str << "{";
		for (int i = 0; i < M.nrows(); i++)
		{
			str << ((i == 0) ? "{ " : " { ");
			for (int j = 0; j < M.ncols(); j++)
				str << M.at(j + i * M.ncols()) << " ";
			str << ((i == M.nrows() - 1) ? "}" : "}\n");
		}
		str << "}\n";
		output << str.str();
		return output;
	}

	template<class T>
	inline typename Shape Matrix<T>::shapeOf(const std::initializer_list<std::initializer_list<T>>& list)
	{
		return Shape{ list.size(), (*list.begin()).size() };
	}

	// Constructors

	template<class T>
	Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& list)
	{
		setShape(shapeOf(list));
		for (auto& row : list)
			for (auto& elem : row)
				matrix_.emplace_back(elem);
	}

	template<class T>
	Matrix<T>::Matrix(const Matrix<T>& other) : matrix_{ other.matrix_ }, shape_{ other.shape_ } {}

	// Shape methods

	template<class T>
	Matrix<T> Matrix<T>::reshape(size_t m, size_t n)
	{
		//if (size() % m || size() % n) throw LinearAlgebraException("Cannot reshape matrix into desired shape.");
		return Matrix{ matrix_, { m, n} };
	}

	template<class T>
	Matrix<T> Matrix<T>::reshape(Shape shape)
	{
		try { reshape(shape.m, shape.n) }
		catch (LinearAlgebraException e) throw e;
	}

	template<class T>
	Matrix<T> Matrix<T>::slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end) const
	{
		size_t m = row_end - row_begin;
		size_t n = column_end - column_begin;
		if (row_end > nrows() || column_end > ncols()) throw LinearAlgebraException("Slice exceeded matrix dimensions.");
		Matrix M = Matrix(m, n);
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				M[i][j] = matrix_[(i + row_begin) * ncols() + j + column_begin];
		return M;
	}

	template<class T>
	void Matrix<T>::insert(Iterator<T> it, const T& val)
	{
		std::ptrdiff_t dist = it - begin();
		std::vector<T>::iterator idx = matrix_.begin() + dist;
		matrix_.insert(idx, val);
	}

	/*
	template<class T>
	void Matrix<T>::insertRows(IteratorRowVector& itIn_beg, IteratorRowVector& itOut_beg,
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
				std::vector<T> v;
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
	template<class T>
	void Matrix<T>::insertColumns(IteratorColumnVector& itIn_beg, IteratorColumnVector& itOut_beg,
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
	template<class T>
	Matrix<T> Matrix<T>::concatenate(const Matrix<T>& other) const
	{
		Matrix<T> mat{ size() + other.size(), { nrows(), ncols() + other.ncols() } };
		for (size_t i = 0; i < mat.nrows(); i++)
		{
			size_t j;
			for (j = 0; j < ncols(); j++)
			{
				mat[i][j] = matrix_[i * ncols() + j];
			}
			for (size_t k = 0; k < other.ncols(); k++, j++)
			{
				mat[i][j] = other.matrix_[i * other.ncols() + k];
			}
		}
		return mat;
	}

	// Iterators and access

	template<class T>
	T& Matrix<T>::at(size_t i)
	{
		if (i > size())
			throw LinearAlgebraException("Index exceeds flattened matrix length");
		return matrix_[i];
	}

	template<class T>
	Iterator<T> Matrix<T>::begin() const { return Iterator<T>{ matrix_.begin()._Ptr }; }

	template<class T>
	Iterator<T> Matrix<T>::end() const { return Iterator<T>{ matrix_.end()._Ptr }; }

template<class T>
Row<T> Matrix<T>::row(size_t i) const
{
	return Row<T>{ begin() + i * ncols(), begin() + (i + 1) * ncols(), getShape() };
}

template<class T>
Column<T> Matrix<T>::col(size_t j) const
{
	return Column<T>
	{
		begin() + j,
		begin() + ncols() * nrows() + j,
		getShape()
	};
}

template<class T>
IteratorRowVector<T> Matrix<T>::beginRow() const { return IteratorRowVector<T>{ begin(), getShape() }; }

template<class T>
IteratorRowVector<T> Matrix<T>::endRow() const { return IteratorRowVector<T>{ begin() + nrows() * ncols(), getShape() }; }

template<class T>
IteratorColumnVector<T> Matrix<T>::beginCol() const { return IteratorColumnVector<T>{ this, 0, getShape() }; }

template<class T>
IteratorColumnVector<T> Matrix<T>::endCol() const { return IteratorColumnVector<T>{ this, ncols(), getShape() }; }

// Algebraic methods
/**/
template<class T> 
Matrix<T> Matrix<T>::inv() const
{
	size_t m = nrows();
	size_t this_n = ncols();
	if (m != this_n) throw LinearAlgebraException("Matrix is not square.");
	Matrix M = Matrix<T>(*this);
	M = M.concatenate(I<T>(m));
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
				T aux = M[i][j];
				M[i][j] = M[pivot][j];
				M[pivot][j] = aux;
			}
		}
		for (size_t k = i + 1; k < m; k++)
		{
			T mki = M[k][i] / M[i][i];
			//M[k][i] = 0;
			for (size_t j = i; j < n; j++)
			{
				M[k][j] -= mki * M[i][j];
			}
		}
	}
	for (int j = m - 1; j >= 0; j--)
	{
		T mjj = 1 / M[j][j];
		for (size_t k = j; k < n; k++)
			M[j][k] *= mjj;
		for (int i = j - 1; i >= 0; i--)
		{
			T mij = M[i][j];
			for (size_t k = j; k < n; k++)
			{
				T mij_mjk = -mij * M[j][k];
				M[i][k] -= mij * M[j][k];
			}
		}
	}
	M = M.slice(0, m, this_n, this_n * 2);
	return M;
}

template<class T>
static T Matrix<T>::determinant_recursion(Matrix M)
{
	size_t m = M.nrows(), n = M.ncols();
	if (m != n)
		throw LinearAlgebraException("Matrix is not square.");
	if (n == 1)
		return M[0];
	if (n == 2)
		return  M[0][0] * M[1][1] - M[0][1] * M[1][0];
	else
	{
		T result = 0;
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

template<class T>
T Matrix<T>::det() const { return determinant_recursion(*this); }

template<class T>
Matrix<T> Matrix<T>::t() const
{
	if (isVector())
	{
		if (isRow())
			return Matrix(matrix_, { ncols(), 1 }); // this constructor creates a column matrix.
		else
			if (isColumn())
			{
				return Matrix(matrix_, { 1, nrows() });
			}
	}
	else
	{
		size_t m = nrows();
		size_t n = ncols();
		Matrix trans{ n, m };
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < m; j++)
				trans[i][j] = matrix_[j * ncols() + i];
		return trans;
	}
	throw LinearAlgebraException("Transpose error: matrix has zero dimensions.");
}
template<class T>
inline T Matrix<T>::norm(T power = 2) const
{
	return sqrt(std::accumulate<Iterator<T>, T>(begin(), end(), 0, [](T& accum, T& next) 
		{ 
			return accum + pow(next, power); 
		}));
}

// Operators

// Unary operator

template<class T>
Matrix<T> Matrix<T>::operator-() const
{
	Matrix<T> mat{ getShape() };
	std::transform(begin(), end(), mat.begin(), std::negate<T>());
	return mat;
}

// Operations with scalars

template<class T>
Matrix<T> Matrix<T>::operator+ (T t) const
{
	Matrix<T> sum{ *this };
	for (T& e : sum) e += t;
	return sum;
}

template<class T>
void Matrix<T>::operator+= (T t) { for (T& e : matrix_) e += t; }

template<class T>
Matrix<T> Matrix<T>::operator- (T t) const 
{
	Matrix<T> sum{ *this };
	for (T& e : sum) e -= t;
	return sum;
}

template<class T>
void Matrix<T>::operator-= (T t) { for (T& e : matrix_) e -= t; }

template<class T>
Matrix<T> Matrix<T>::operator* (T t) const
{
	Matrix mul{ getShape() };
	for (T& e : mul) e *= t;
	return mul;
}

template<class T>
void Matrix<T>::operator*= (T t) { for (T& e : matrix_) e *= t; }

template<class T>
Matrix<T> Matrix<T>::operator/ (T t) const
{
	Matrix mul{ *this }; 
	for (T& e : mul) e /= t;
	return mul;
}

template<class T>
void Matrix<T>::operator/= (T t) { for (T& e : matrix_) e /= t; }

// Operations with Matrix

template<class T>
inline Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const
{
	if (nrows() != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
	if (ncols() != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
	Matrix<T> sum{ getShape() };
	std::transform(begin(), end(), other.begin(), sum.begin(), std::plus<T>());
	return sum;
}

template<class T>
inline void Matrix<T>::operator+=(const Matrix<T>& other)
{
	std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
}

template<class T>
inline Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const
{
	if (nrows() != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
	if (ncols() != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
	Matrix<T> sum{ getShape() };
	std::transform(begin(), end(), other.begin(), sum.begin(), std::minus<T>());
	return sum;
}

template<class T>
inline void Matrix<T>::operator-=(const Matrix<T>& other)
{
	std::transform(begin(), end(), other.begin(), begin(), std::negate<T>());
}

template<class T>
inline Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const
{
	size_t m = nrows();
	size_t n = ncols();
	size_t q = other.ncols();
	if (other.isScalar())
		return (*this) * other.matrix_[0];
	else if (isScalar())
		return other * matrix_[0];
	else if (n != other.nrows()) throw LinearAlgebraException("Multiplication error: matrices do not have appropriate dimensions.");
	else
	{
		//Matrix mul(m, q);
		Matrix mul(0, { m, q });
		std::for_each(beginRow(), endRow(), [&](Row<T>& row)
			{
				std::for_each(other.beginCol(), other.endCol(), [&](Column<T>& col)
					{
						mul.matrix_.push_back(
							std::inner_product(row.begin(), row.end(), col.begin(), 0)
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

template<class T>
inline void Matrix<T>::operator*=(const Matrix<T>& other)
{
	size_t m = nrows();
	if (m != ncols())
		throw LinearAlgebraException("Multiplication error: assigning result of multiplication between non-square matrices to self matrix.");
	if (m != other.nrows() || m != other.ncols())
		throw LinearAlgebraException("Multiplication error: assining result of non-square matrix multiplication to left-hand-side matrix.");
	(*this) = (*this) * other;
}




/**/
} // namespace::alg