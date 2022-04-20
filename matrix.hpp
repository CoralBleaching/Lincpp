#pragma once

#include<vector>
#include<algorithm>
#include<numeric>
#include<sstream>
#include<ostream>
#include<cmath>

namespace alg {

template<class T> class Matrix;
template<class T> Matrix<T> I(size_t);
class LinearAlgebraException;
template<typename Container> std::string toString(Container con);

template<class T>
class Matrix
{

	// DECLARATIONS
private:

	static T determinant_recursion(Matrix<T>);

protected:
	std::vector<T> matrix_;
	struct Shape { size_t m, n; } shape_;

	static inline Shape shapeOf(const std::initializer_list<std::initializer_list<T>>& );

public: // PRELOAD ALL ROWS AND COLUMNS FOR SPEED EFFICIENCY?
	Matrix(size_t k = 0, Shape shape = { 0, 0 }) :
		matrix_{ std::vector<T>(k, 0) },
		shape_{ shape } 
	{}
	Matrix(size_t m, size_t n) : matrix_{ std::vector<T>(m * n, 0) }, shape_{ m, n } {}
	Matrix(Shape shape) : matrix_{ std::vector<T>(m * n, 0) }, shape_{ shape } {}
	Matrix(std::vector<T> matrix, Shape shape) : matrix_{ matrix }, shape_{ shape } {}
	Matrix(std::vector<T> matrix, size_t m, size_t n) : matrix_{ matrix }, shape_{ m, n } {}
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
	inline bool isScalar() const { return nrows() == ncols() == 1; }
	inline bool isRow() const { return ncols() == 1; }
	inline bool isColumn() const { return nrows() == 1; }
	inline bool isVector() const { return isRow() || isColumn(); }
	Matrix<T> slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end) const;
	inline void clear() { matrix_.clear(); shape_ = { 0, 0 }; }

	class Iterator;
	class IteratorColumn;
	class Row;
	class Column;
	class IteratorRowVector;
	class IteratorColumnVector;

	void insert(Iterator it, const T& val);
	//void insertRows(IteratorRowVector& itIn_beg, IteratorRowVector& itOut_beg, IteratorRowVector& itOut_end);
	//void insertColumns(IteratorColumnVector& itIn_beg, IteratorColumnVector& itOut_beg, IteratorColumnVector& itOut_end);
	Matrix<T> concatenate(const Matrix<T>&);
	//Matrix<T> concatenateBelow(const Matrix<T>&);
	//virtual Matrix<T> slice(size_t index_begin, size_t index_end) const;
	//virtual Matrix<T> slice(const arrayIterator<T>& it_b, const arrayIterator<T>& it_e) const;

	// Iterators and access

	inline Iterator begin() const;
	inline Iterator end() const;
	inline IteratorRowVector beginRow() const;
	inline IteratorRowVector endRow() const;
	inline IteratorColumnVector beginCol() const;
	inline IteratorColumnVector endCol() const;

	T& at(size_t i);
	Row operator[](size_t i) const;
	Column col(size_t j) const;

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
	T norm() const;
	friend T norm(Matrix<T> M) { return M.norm(); }
	// Operators

	// Unary operator

	//Matrix<T> operator- () const;

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

	virtual Matrix<T> operator+ (const Matrix<T>&) const;
	virtual void operator+= (const Matrix<T>&);
	virtual Matrix<T> operator- (const Matrix<T>&) const;
	virtual void operator-= (const Matrix<T>&);

	Matrix<T> operator* (const Matrix<T>&) const;
	void operator*= (const Matrix<T>&);

/**/
	// DEFINITIONS


	class Iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using reference = T&;
		Iterator(pointer ptr = nullptr) : mptr{ ptr } {}
		Iterator(const Iterator& rawIterator) = default;
		Iterator& operator=(const Iterator& rawIterator) = default;
		reference operator*() const { return *mptr; }
		Iterator& operator++() { mptr++; return *this; }
		Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
		Iterator& operator--() { mptr--; return *this; }
		Iterator operator--(int) { Iterator tmp = *this; --(*this); return tmp; }
		Iterator operator+(difference_type movement)
		{
			auto oldptr = mptr; mptr += movement; auto tmp{ *this }; mptr = oldptr;
			return tmp;
		}
		Iterator operator-(difference_type movement)
		{
			auto oldptr = mptr; mptr -= movement; auto tmp{ *this }; mptr = oldptr;
			return tmp;
		}
		void operator+=(difference_type movement) { mptr += movement; }
		void operator-=(difference_type movement) { mptr -= movement; }
		friend difference_type operator-(const Iterator& it1, const Iterator& it2)
		{
			return std::distance(it2.getPtr(), it1.getPtr());
		}
		friend bool operator==(const Iterator& a, const Iterator& b) { return a.mptr == b.mptr; }
		friend bool operator!=(const Iterator& a, const Iterator& b) { return a.mptr != b.mptr; }

		inline T* getPtr() const { return mptr; }

	private:
		T* mptr;
	};

	class IteratorColumn
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using reference = T&;

		IteratorColumn(pointer ptr, size_t ncols) : ncols_{ ncols } { mptr = ptr; }
		IteratorColumn(Iterator it, size_t ncols) : ncols_{ ncols } { mptr = it.getPtr(); }
		IteratorColumn(const IteratorColumn& it) : ncols_{ it.ncols_ } { mptr = it.getPtr(); }
		reference operator*() const { return *mptr; }
		IteratorColumn& operator++() { mptr += ncols_; return *this; }
		IteratorColumn operator++(int) { IteratorColumn tmp = *this; ++(*this); return tmp; }
		IteratorColumn& operator--() { mptr -= ncols_; return *this; }
		IteratorColumn operator--(int) { IteratorColumn tmp = *this; --(*this); return tmp; }
		IteratorColumn operator+(difference_type movement)
		{
			auto oldptr = mptr; mptr += movement * ncols_; auto tmp{ *this }; mptr = oldptr;
			return tmp;
		}
		IteratorColumn operator-(difference_type movement)
		{
			auto oldptr = mptr; mptr -= movement * ncols_; auto tmp{ *this }; mptr = oldptr;
			return tmp;
		}
		void operator+=(difference_type movement) { mptr += movement * ncols_; }
		void operator-=(difference_type movement) { mptr -= movement * ncols_; }
		friend difference_type operator-(const IteratorColumn& it1, const IteratorColumn& it2) // relative distance
		{
			return std::distance(it2.getPtr(), it1.getPtr()) / it1.ncols();
		}
		friend difference_type abs_dist(const IteratorColumn& it1, const IteratorColumn& it2) // absolute distance
		{
			return std::distance(it2.getPtr(), it1.getPtr())();
		}

		friend bool operator==(const IteratorColumn& a, const IteratorColumn& b) { return a.mptr == b.mptr; }
		friend bool operator!=(const IteratorColumn& a, const IteratorColumn& b) { return a.mptr != b.mptr; }

		inline T* getPtr() const { return mptr; }
		inline size_t ncols() const { return ncols_; }

	private:
		T* mptr;
		size_t ncols_;
	};

public:

	class Row
	{
	public:
		Row(Iterator begin = Iterator{ nullptr },
			Iterator end = Iterator{ nullptr },
			Shape shape = { 0,0 }) :
			begin_{ begin }, end_{ end }, shape_{ shape } {} // automatically generate end_?

		Iterator begin() const { return begin_; }
		Iterator end() const { return end_; }

		inline Shape getShape() { return shape_; }
		inline size_t ncols() { return shape_.n; }
		inline size_t size() { return shape_.n; }

		T& operator[](size_t i) const { return *(begin_.getPtr() + i); } //THROW
		friend std::ostream& operator<< (std::ostream& output, Row r) {	output << toString(r); return output; }

	protected:
		Iterator begin_, end_;
		Shape shape_;
	};

	class Column
	{
	public:
		Column(Iterator begin = Iterator{ nullptr },
			Iterator end = Iterator{ nullptr },
			Shape shape = { 0,0 }) :
			begin_{ IteratorColumn{begin, shape.n} }, end_{ IteratorColumn{end, shape.n} }, shape_{ shape } {}

		IteratorColumn begin() const { return begin_; }
		IteratorColumn end() const { return end_; }

		inline Shape getShape() { return shape_; }
		inline size_t nrows() { return shape_.m; }
		inline size_t size() { return shape_.m; }

		T& operator[](size_t i) { return *(begin_ + i); } //THROW
		friend std::ostream& operator<< (std::ostream& output, Column c) { output << toString(c); return output; }

	private:
		IteratorColumn begin_, end_;
		Shape shape_;
	};

public:

	class IteratorRowVector
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		//using value_type = Row;
		//using pointer = Iterator;
		//using reference = Row&;

		IteratorRowVector(Iterator it = Iterator{ nullptr }, Shape shape = { 0,0 }) : it_{ it }, shape_{ shape } {}
		IteratorRowVector(const IteratorRowVector& rawIterator) = default;
		IteratorRowVector& operator=(const IteratorRowVector& rawIterator) = default;

		Row& operator*() { Row r{ it_, it_ + ncols(), shape_ }; return r; }

		IteratorRowVector& operator++() { it_ += ncols(); return *this; }
		IteratorRowVector operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		IteratorRowVector& operator--() { it_ -= ncols(); return *this; }
		IteratorRowVector operator--(int) { auto tmp = *this; --(*this); return tmp; }
		IteratorRowVector operator+(difference_type movement)
		{
			auto olditr = it_; it_ += movement * ncols(); auto tmp{ *this }; it_ = olditr;
			return tmp;
		}
		IteratorRowVector operator-(difference_type movement)
		{
			auto olditr = it_; it_ -= movement * ncols(); auto tmp{ *this }; it_ = olditr;
			return tmp;
		}
		void operator+=(difference_type movement) { it_ += movement * ncols(); }
		void operator-=(difference_type movement) { it_ -= movement * ncols(); }
		friend bool operator==(const IteratorRowVector& a, const IteratorRowVector& b) { return a.it_ == b.it_; }
		friend bool operator!=(const IteratorRowVector& a, const IteratorRowVector& b) { return a.it_ != b.it_; }

		Iterator getIt() { return it_; }
		size_t nrows() { return shape_.m; }
		size_t ncols() { return shape_.n; }

	private:
		Iterator it_;
		Shape shape_;
	};

	class IteratorColumnVector
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;

		IteratorColumnVector(const Matrix<T>* ptr = nullptr, size_t index = 0, Shape shape = { 0,0 }) : mptr{ ptr }, index_{ index }, shape_{ shape } {}

		Column& operator*()
		{
			Column c = mptr->col(index_);
			return c;
		}
		IteratorColumnVector& operator++() { index_ += 1; return *this; }
		IteratorColumnVector operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		IteratorColumnVector& operator--() { index_ -= 1; return *this; }
		IteratorColumnVector operator--(int) { auto tmp = *this; --(*this); return tmp; }
		IteratorColumnVector operator+(difference_type movement)
		{
			auto oldidx = index_; index_ += movement; auto tmp{ *this }; index_ = oldidx;
			return tmp;
		}
		IteratorColumnVector operator-(difference_type movement)
		{
			auto oldidx = index_; index_ -= movement; auto tmp{ *this }; index_ = oldidx;
			return tmp;
		}
		void operator+=(difference_type movement) { index_ += movement; }
		void operator-=(difference_type movement) { index_ -= movement; }

		friend bool operator==(const IteratorColumnVector& a, const IteratorColumnVector& b)
		{
			return a.index_ == b.index_ && a.mptr == b.mptr;
		}
		friend bool operator!=(const IteratorColumnVector& a, const IteratorColumnVector& b) { return !(a == b); }

		Matrix<T>* getMatrixPtr() { return mptr; }
		size_t getIndex() { return index_; }
		size_t nrows() { return shape_.m; }
		size_t ncols() { return shape_.n; }

	private:
		const Matrix<T>* mptr;
		size_t index_;
		Shape shape_;
	};

};

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
inline typename Matrix<T>::Shape Matrix<T>::shapeOf(const std::initializer_list<std::initializer_list<T>>& list)
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
void Matrix<T>::insert(Iterator it, const T& val)
{
	std::ptrdiff_t dist = it - begin();
	std::vector<T>::iterator idx = matrix_.begin() + dist;
	matrix_.insert(idx, val);
}

/**
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
Matrix<T> Matrix<T>::concatenate(const Matrix<T>& other)
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
typename Matrix<T>::Iterator Matrix<T>::begin() const { return Iterator{ matrix_.begin()._Ptr }; }

template<class T>
typename Matrix<T>::Iterator Matrix<T>::end() const { return Iterator{ matrix_.end()._Ptr }; }

template<class T>
typename Matrix<T>::Row Matrix<T>::operator[](size_t i) const
{
	return Row{ begin() + i * ncols(), begin() + (i + 1) * ncols(), getShape() };
}

template<class T>
typename Matrix<T>::Column Matrix<T>::col(size_t j) const
{
	return Column
	{
		begin() + j,
		begin() + ncols() * nrows() + j,
		getShape()
	};
}

template<class T>
typename Matrix<T>::IteratorRowVector Matrix<T>::beginRow() const { return IteratorRowVector{ begin(), getShape() }; }

template<class T>
typename Matrix<T>::IteratorRowVector Matrix<T>::endRow() const { return IteratorRowVector{ end(), getShape() }; }

template<class T>
typename Matrix<T>::IteratorColumnVector Matrix<T>::beginCol() const { return IteratorColumnVector{ this, 0, getShape() }; }

template<class T>
typename Matrix<T>::IteratorColumnVector Matrix<T>::endCol() const { return IteratorColumnVector{ this, ncols(), getShape() }; }

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
		return M[0] * M[1 * ncols() + 1] - M[0 * ncols() + 1] * M[1 * ncols() + 0];
	else
	{
		T result = 0;
		for (size_t i = 0; i < n; i++)
		{
			Matrix left_submatrix = M.slice(1, n, 0, i);
			Matrix right_submatrix = M.slice(1, n, i + 1, n);
			Matrix submatrix = left_submatrix.concatenate(right_submatrix);
			result += std::pow(-1, i) * M[i] * determinant_recursion(submatrix);
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
			return Matrix(matrix_, { nrows(), 1 }); // this constructor creates a column matrix.
		else
			if (isColumn())
			{
				//size_t n = M.n_rows();
				//Matrix matrix = Matrix<T>(1, n);
				//for (size_t j = 0; j < n; j++) matrix[0][j] = M[j][0];
				return Matrix(matrix_, { 1, ncols() });
			}
	}
	else
	{
		Matrix matrix(*this);
		size_t m = nrows();
		size_t n = ncols();
		if (m != n)
		{
			std::ostringstream error_msg("Transpose error: non-square matrix. Dimensions: ");
			error_msg << m << "x" << n << std::endl;
			throw LinearAlgebraException(error_msg.str().data());
		}
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < i; j++)
			{
				T aux = matrix[i * ncols() + j];
				matrix[i * ncols() + j] = matrix[j * ncols() + i];
				matrix[j * ncols() + i] = aux;
			}
		return matrix;
	}
	throw LinearAlgebraException("Transpose error: matrix has zero dimensions.");
}
template<class T>
inline T Matrix<T>::norm() const
{
	return sqrt(std::accumulate(begin(), end(), 0, [](T& accum, T& next) { return accum + pow(next, 2) });
}

// Operations with scalars

template<class T>
Matrix<T> Matrix<T>::operator+ (T t) const
{
	Matrix<T> sum = Matrix<T>(*this);
	size_t m = n_rows();
	size_t n = n_columns();
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			sum[i][j] += t;
	return sum;
}

template<class T>
void Matrix<T>::operator+= (T t)
{
	*this = *this + t;
}

template<class T>
Matrix<T> Matrix<T>::operator- (T t) const
{
	Matrix<T> sum = Matrix<T>(*this);
	size_t m = n_rows();
	size_t n = n_columns();
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			sum[i][j] += t;
	return sum;
}

template<class T>
void Matrix<T>::operator-= (T t)
{
	*this = *this - t;
}

template<class T>
Matrix<T> Matrix<T>::operator* (T t) const
{
	size_t m = nrows();
	size_t n = ncols();
	Matrix mul{ *this };
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			mul[i * ncols() + j] *= t;
	return mul;
}

template<class T>
void Matrix<T>::operator*= (T t)
{
	*this = *this * t;
}

template<class T>
Matrix<T> Matrix<T>::operator/ (T t) const
{
	size_t m = nrows();
	size_t n = ncols();
	Matrix mul{ *this };
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			mul[i * ncols() + j] /= t;
	return mul;
}

template<class T>
void Matrix<T>::operator/= (T t)
{
	*this = *this / t;
}
// Operations with Matrix

template<class T>
inline Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const
{
	size_t m = nrows();
	size_t n = ncols();
	if (m != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
	if (n != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
	Matrix<T> sum = Matrix<T>(m, n);
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			sum[i][j] = matrix_[i * ncols() + j] + other.matrix_[i * ncols() + j];
	return sum;
}

template<class T>
inline void Matrix<T>::operator+=(const Matrix<T>& other)
{
	*this = *this + other;
}

template<class T>
inline Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const
{
	size_t m = nrows();
	size_t n = ncols();
	if (m != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
	if (n != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
	Matrix<T> sum = Matrix<T>(m, n);
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			sum[i][j] = matrix_[i * ncols() + j] - other.matrix_[i * ncols() + j];
	return sum;
}

template<class T>
inline void Matrix<T>::operator-=(const Matrix<T>& other)
{
	*this = *this - other;
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
		Matrix mul{ m, q };
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < q; j++)
				for (size_t k = 0; k < n; k++)
					mul[i][j] += matrix_[i * ncols() + k] * other.matrix_[k * ncols() + j];
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