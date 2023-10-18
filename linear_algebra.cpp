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

	Vector::Vector(size_type k, value_type val) : data_{ std::vector<value_type>(k, val) } {}
	Vector::Vector(const std::vector<value_type> v) : data_{ v } {}
	Vector::Vector(const std::initializer_list<value_type> l) : data_{ l } {}
	Vector::Vector(const_iterator first, const_iterator last) : data_{ std::vector<value_type>(first, last) } {}
	Vector::Vector(Matrix::const_column_iterator first, Matrix::const_column_iterator last) : data_{ std::vector<value_type>(first, last) } {}
	Vector::Vector(const Matrix& M) : data_{ M.getInternalStdVector() } {}

	// Shape methods

	Vector::size_type Vector::size() const { return data_.size(); }
	bool              Vector::isScalar() const { return size() == 1; }
	void              Vector::clear() { data_.clear(); }

	Matrix Vector::row() const { return Matrix{ data_, 1, size() }; }
	Matrix row(Vector v) { return v.row(); }
	Matrix Vector::col() const { return Matrix{ data_, size(), 1 }; }
	Matrix col(Vector v) { return v.col(); }
	Matrix Vector::t() const { return row(); }
	Matrix t(Vector v) { return v.t(); }

	// Iterators and access

	std::vector<Vector::value_type>& Vector::getInternalStdVector() { return data_; }

	void Vector::push_back(value_type val) { data_.push_back(val); }

	void Vector::insert(iterator it, value_type val)
	{
		std::ptrdiff_t dist = it - begin();
		std::vector<value_type>::iterator pos = data_.begin() + dist;
		data_.insert(pos, val);
	}

	void Vector::insert(iterator it, const Vector& v)
	{
		std::ptrdiff_t dist = it - begin();
		std::vector<value_type>::iterator pos = data_.begin() + dist;
		std::copy(v.begin(), v.end(), std::inserter(data_, pos));
	}

	Vector Vector::concatenate(const Vector& v) const
	{
		Vector res(*this);
		res.insert(res.end(), v);
		return res;
	}

	Vector Vector::slice(size_type start, size_type finish)
	{
		size_type n = finish - start;
		if (n > size()) throw LinearAlgebraException("Slice exceeds Vector length.");
		Vector v(n);
		std::copy_n(begin() + start, n, v.begin());
		return v;
	}

	Vector::iterator Vector::begin() { return iterator{ data_.data() }; }
	Vector::iterator Vector::end() { return iterator{ data_.data() + data_.size() }; }
	Vector::const_iterator Vector::begin() const { return const_iterator(data_.data()); }
	Vector::const_iterator Vector::end() const { return const_iterator(data_.data() + data_.size()); }

	Vector::reference Vector::operator[](size_type i) 
	{ 
		if (i >= size()) throw LinearAlgebraException("Invalid access: index exceeds Vector length.");
		return data_[i]; 
	}
	Vector::reference Vector::at(size_type i) 
	{ 
		if (i >= size()) throw LinearAlgebraException("Invalid access: index exceeds Vector length.");
		return data_.at(i); 
	}

	std::string Vector::to_string() const
	{
		std::ostringstream oss;
		oss << *this;
		return oss.str();
	}

	// Algebraic methods

	Vector::value_type Vector::norm(value_type power) const
	{
		return sqrt(std::accumulate(begin(), end(), static_cast<value_type>(0), [=](value_type accum, value_type next)
			{
				return accum + pow(next, power);
			}));
	}

	Vector::value_type norm(Vector v, Vector::value_type val) { return v.norm(val); }

	// Operators

	// Unary operator

	Vector Vector::operator-() const
	{
		Vector v(size());
		std::transform(begin(), end(), v.begin(), std::negate<value_type>());
		return v;
	}

	// Operations with scalar

	Vector Vector::operator+(value_type t) const
	{
		Vector v{ *this };
		for (reference e : v) e += t;
		return v;
	}

	void Vector::operator+=(value_type t) { for (reference e : data_) e += t; }

	Vector operator+(Vector::value_type t, Vector v) { return v + t; }

	Vector Vector::operator-(value_type t) const
	{
		Vector v{ *this };
		for (reference e : v) e -= t;
		return v;
	}

	void Vector::operator-=(value_type t) { for (reference e : data_) e -= t; }

	Vector operator-(Vector::value_type t, Vector v)
	{
		for (Vector::reference e : v) e = t - e;
		return v;
	}

	Vector Vector::operator*(value_type t) const
	{
		Vector v{ *this };
		for (reference e : v) e *= t;
		return v;
	}

	void Vector::operator*=(value_type t) { for (reference e : data_) e *= t; }

	Vector operator*(Vector::value_type t, Vector v) { return v * t; }

	Vector Vector::operator/(value_type t) const
	{
		Vector v{ *this };
		for (reference e : v) e /= t;
		return v;
	}
	void Vector::operator/=(value_type t) { for (reference e : data_) e /= t; }

	Vector operator/(Vector::value_type t, Vector v) { for (Vector::reference e : v) e = t / e; return v; }

	// Operations with Vector

	Vector::value_type Vector::operator*(Vector v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}

	Vector Vector::operator+(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}

	void Vector::operator+=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}

	Vector Vector::operator-(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}

	void Vector::operator-=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}


	// Operations with Matrix


	Vector Vector::operator+(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}

	void Vector::operator+=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}

	Vector Vector::operator-(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}

	void Vector::operator-=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	Vector Vector::operator*(Matrix M) const
	{
		if (size() != M.nrows()) throw LinearAlgebraException("Can't perform dot product: Vector and Matrix have different number of columns.");
		Vector v;
		std::for_each(M.beginCol(), M.endCol(), [&](Matrix::Column col) // make const iterator to allow & argument
			{
				v.push_back(std::inner_product(col.begin(), col.end(), begin(), 0.));
			});
		return v;
	}

	// Operations with Row, Column

	Vector::value_type Vector::operator*(const Matrix::Row& r) const { return r * (*this); }

	Vector Vector::operator+(const Matrix::Row& r) const { return r + *this; }

	Vector Vector::operator-(const Matrix::Row& r) const
	{
		if (r.size() != size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v;
		std::transform(begin(), end(), r.begin(), std::back_inserter(v.data_), std::minus<value_type>());
		return v;
	}

	void Vector::operator+=(const Matrix::Row& r)
			{
		std::transform(begin(), end(), r.begin(), begin(), std::plus<value_type>());
	}

	void Vector::operator-=(const Matrix::Row& r)
	{
		std::transform(begin(), end(), r.begin(), begin(), std::minus<value_type>());
	}

	Vector::value_type Vector::operator*(const Matrix::Column& r) const { return r * (*this); }

	Vector Vector::operator+(const Matrix::Column& r) const { return r + *this; }

	Vector Vector::operator-(const Matrix::Column& r) const
	{
		if (r.size() != size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v;
		std::transform(begin(), end(), r.begin(), std::back_inserter(v.data_), std::minus<value_type>());
		return v;
	}

	void Vector::operator+=(const Matrix::Column& r)
	{
		std::transform(begin(), end(), r.begin(), begin(), std::plus<value_type>());
	}

	void Vector::operator-=(const Matrix::Column& r)
	{
		std::transform(begin(), end(), r.begin(), begin(), std::minus<value_type>());
	}


	// Other operations

	std::ostream& operator<<(std::ostream& ostream, const Vector& vector)
	{
		ostream << "{ ";
		for (auto& elem : vector) ostream << elem << " ";
		ostream << "}";
		return ostream;
	}

	void Vector::operator=(const Matrix& m)
	{
		if (!m.isVector()) throw LinearAlgebraException("Can't assign to Vector: Matrix isn't row or column.");
		if (size() != m.length()) throw LinearAlgebraException("Can't assign Matrix to Vector of different length.");
		data_.assign(m.begin(), m.end());
	}

	void Vector::operator=(const Matrix::Row& v)
	{
		if (size() != v.size()) throw LinearAlgebraException("Can't assign to Vector of different length.");
		data_.assign(v.begin(), v.end());
	}

	void Vector::operator=(const Matrix::Column& v)
	{
		if (size() != v.size()) throw LinearAlgebraException("Can't assign to Vector of different length.");
		data_.assign(v.begin(), v.end());
	}


	// MATRIX


	// Utils

	Matrix I(Matrix::size_type m)
	{
		Matrix e{ m, m };
		for (Matrix::size_type i = 0; i < m; i++) e[i][i] = static_cast<Matrix::value_type>(1);
		return e;
	}

	Vector arange(Matrix::value_type end, Matrix::value_type start, Matrix::value_type step)
	{
		if ((start > end && step > 0.) || (start < end && step < 0.) || step == 0.)
			throw LinearAlgebraException("Invalid arguments to function range().");
		Vector rangelist(static_cast<Matrix::size_type>(ceil((end - start)) / step), start);
		Matrix::value_type loopstep = step;
		std::transform(rangelist.begin() + 1, rangelist.end(), rangelist.begin() + 1, [&](Matrix::reference e) {
			Matrix::value_type tmp = e + loopstep; loopstep += step; return tmp;
			});
		return rangelist;
	}

	Shape Matrix::shapeOf(const std::initializer_list<std::initializer_list<Matrix::value_type>>& list)
	{
		if (list.size() < 1) return { 0, 0 };
		return { list.size(), (*list.begin()).size() };
	}


	// Constructors

	alg::Matrix::Matrix(alg::Shape shape, value_type val) : data_{ std::vector<value_type>(shape.m * shape.n, val) }, shape_{ shape } {}

	alg::Matrix::Matrix(size_type m, size_type n, value_type val) : data_{ std::vector<value_type>(m * n, val) }, shape_{ m, n } {}

	alg::Matrix::Matrix(std::vector<value_type> data, Shape shape) : data_{ data }, shape_{ shape } 
	{
		if (data.size() != shape.m * shape.n) {
			std::ostringstream msg;
			msg << "Can't conform std::vector of size "
				<< data.size() << " to alg::Matrix of "
				<< "shape " << shape.m << "x" << shape.n << ".";
			throw LinearAlgebraException(msg.str());
		}
	}

	alg::Matrix::Matrix(std::vector<value_type> data, size_type m, size_type n) : data_{ data }, shape_{ m, n } 
	{
		if (data.size() != m * n) {
			std::ostringstream msg;
			msg << "Can't conform std::vector of size "
				<< data.size() << " to alg::Matrix of "
				<< "shape " << m << "x" << n << ".";
			throw LinearAlgebraException(msg.str());
		}
	}

	alg::Matrix::Matrix(std::vector<value_type> data) : data_{ data }, shape_{ data.size(), 1 } {}

	alg::Matrix::Matrix(Vector vector) : data_{ vector.getInternalStdVector() }, shape_{ vector.size(), 1 } {}

	alg::Matrix::Matrix(alg::Vector vector, Shape shape) : data_{ vector.getInternalStdVector() }, shape_{ shape } 
	{
		if (vector.size() != shape.m * shape.n) {
			std::ostringstream msg;
			msg << "Can't conform alg::Vector of size "
				<< vector.size() << " to alg::Matrix of "
				<< "shape " << shape.m << "x" << shape.n << ".";
			throw LinearAlgebraException(msg.str());
		}
	}

	alg::Matrix::Matrix(alg::Vector vector, size_type m, size_type n) : data_{ vector.getInternalStdVector() }, shape_{ m, n } 
	{
		if (vector.size() != m * n) {
			std::ostringstream msg;
			msg << "Can't conform alg::Vector of size "
				<< vector.size() << " to alg::Matrix of "
				<< "shape " << m << "x" << n << ".";
			throw LinearAlgebraException(msg.str());
		}
	}

	alg::Matrix::Matrix(const std::initializer_list<std::initializer_list<value_type>>& list)
	{
		setShape(shapeOf(list));
		for (auto& row : list)
			for (auto& elem : row)
				data_.emplace_back(elem);
	}

	// Shape methods


	Matrix::size_type Matrix::length() const { return data_.size(); }
	Shape Matrix::getShape() const { return shape_; }
	Matrix::size_type Matrix::nrows() const { return shape_.m; }
	Matrix::size_type Matrix::ncols() const { return shape_.n; }

	void Matrix::setShape(size_type m, size_type n) 
	{ 
		setShape({m, n}); 
	}

	void Matrix::setShape(Shape shape) 
	{ 
		if (length() != shape.m * shape.n) 
		{
			std::ostringstream msg;
			msg << "Can't conform " << nrows() << "x"
				<< ncols() << " alg::Matrix to "
				<< shape.m << "x" << shape.n << ".";
			throw LinearAlgebraException(msg.str());
		}
		shape_ = shape; 
	}

	Matrix Matrix::reshape(size_type m, size_type n) const 
	{
		if (length() % m || length() % n) 
		{
			std::ostringstream msg;
			msg << "Can't conform " << nrows() << "x"
				<< ncols() << " alg::Matrix to "
				<< m << "x" << n << ".";
			throw LinearAlgebraException(msg.str());
		}
		return Matrix{ data_, { m, n } };
	};

	Matrix Matrix::reshape(Shape shape) const
	{
		return reshape(shape.m, shape.n);
	}

	bool Matrix::isScalar() const { return nrows() == 1 && ncols() == 1; }
	bool Matrix::isRow() const { return nrows() == 1; }
	bool Matrix::isColumn() const { return ncols() == 1; }
	bool Matrix::isVector() const { return isRow() || isColumn(); }
	void Matrix::clear() { data_.clear(); shape_ = { 0, 0 }; }

	// Insertions and slice

	void Matrix::insertRow(const Vector& v)
	{
		insertRow(v.begin(), v.end());
	}

	void Matrix::insertRow(const Vector& v, size_type rowIndex)
	{
		insertRow(v.begin(), v.end(), rowIndex);
	}

	void Matrix::insertRow(const Matrix::Row& v)
	{
		insertRow(v.begin(), v.end());
	}

	void Matrix::insertRow(const Matrix::Row& v, size_type rowIndex)
	{
		insertRow(v.begin(), v.end(), rowIndex);
	}

	void Matrix::insertRow(const Matrix::Column& v)
	{
		insertRow(v.begin(), v.end());
	}

	void Matrix::insertRow(const Matrix::Column& v, size_type rowIndex)
	{
		insertRow(v.begin(), v.end(), rowIndex);
	}

	void Matrix::insertRow(const Matrix& v)
	{
		insertRow(v.begin(), v.end());
	}

	void Matrix::insertRow(const Matrix& v, size_type rowIndex)
	{
		insertRow(v.begin(), v.end(), rowIndex);
	}

	void Matrix::insertRow(const_iterator begin, const_iterator end)
	{
		insertRow(begin, end, shape_.m);
	}

	void Matrix::insertRow(const_column_iterator begin, const_column_iterator end)
	{
		insertRow(begin, end, shape_.m);
	}
	
	void Matrix::insertRow(const_iterator begin, const_iterator end, size_type rowIndex)
	{
		difference_type length = end - begin;
		if (nrows() > 0)
			if (length != ncols()) 
				throw LinearAlgebraException("Can't insert row of incompatible length.");
		else
			setShape(0, length);
		std::vector<value_type>::iterator idx = data_.begin() + rowIndex * length;
		data_.insert(idx, begin, end);
		setShape(nrows() + 1, ncols());
	}
	
	void Matrix::insertRow(const_column_iterator begin, const_column_iterator end, size_type rowIndex)
	{
		difference_type length = end - begin;
		if (nrows() > 0)
			if (length != ncols()) 
				throw LinearAlgebraException("Can't insert row of incompatible length.");
		else
			setShape(0, length);
		std::vector<value_type>::iterator idx = data_.begin() + rowIndex * length;
		data_.insert(idx, begin, end);
		setShape(nrows() + 1, ncols());
	}

	void Matrix::insertColumn(const Vector& v)
	{
		insertColumn(v.begin(), v.end());
	}

	void Matrix::insertColumn(const Vector& v, size_type columnIndex)
	{
		insertColumn(v.begin(), v.end(), columnIndex);
	}

	void Matrix::insertColumn(const Matrix::Row& v)
	{
		insertColumn(v.begin(), v.end());
	}

	void Matrix::insertColumn(const Matrix::Row& v, size_type columnIndex)
	{
		insertColumn(v.begin(), v.end(), columnIndex);
	}

	void Matrix::insertColumn(const Matrix::Column& v)
	{
		insertColumn(v.begin(), v.end());
	}

	void Matrix::insertColumn(const Matrix::Column& v, size_type columnIndex)
	{
		insertColumn(v.begin(), v.end(), columnIndex);
	}

	void Matrix::insertColumn(const_column_iterator begin, const_column_iterator end)
	{
		insertColumn(begin, end, shape_.n);
	}

	void Matrix::insertColumn(const_iterator begin, const_iterator end)
	{
		insertColumn(begin, end, shape_.n);
	}

	void Matrix::insertColumn(const_column_iterator begin, const_column_iterator end, size_type columnIndex)
	{
		difference_type length = end - begin;
		if (ncols() > 0)
			if (length != nrows()) 
				throw LinearAlgebraException("Can't insert column of incompatible length.");
		else
			setShape(length, 0);
		difference_type count = 0, offset = 0;
		std::for_each(begin, end, [&](value_type e) {
			std::vector<value_type>::iterator idx = data_.begin() + columnIndex + shape_.n * count++ + offset++;
			data_.insert(idx, e);
		});
		shape_.n++;
	}

	void Matrix::insertColumn(const_iterator begin, const_iterator end, size_type columnIndex)
	{
		difference_type length = end - begin;
		if (ncols() > 0)
			if (length != nrows()) 
				throw LinearAlgebraException("Can't insert column of incompatible length.");
		else
			setShape(length, 0);
		difference_type count = 0, offset = 0;
		std::for_each(begin, end, [&](value_type e) {
			std::vector<value_type>::iterator idx = data_.begin() + columnIndex + shape_.n * count++ + offset++;
			data_.insert(idx, e);
		});
		shape_.n++;
	}

	Matrix Matrix::slice(size_type row_begin, size_type row_end, size_type column_begin, size_type column_end) const
	{
		size_type m = row_end - row_begin;
		size_type n = column_end - column_begin;
		if (row_end > nrows() || column_end > ncols()) throw LinearAlgebraException("Slice exceeded matrix dimensions.");
		Matrix M = Matrix(m, n);
		for (size_type i = 0; i < m; i++)
			for (size_type j = 0; j < n; j++)
				M(i, j) = data_[(i + row_begin) * ncols() + j + column_begin];
		return M;
	}


	Matrix Matrix::concatenate(const Matrix& other) const
	{
		Matrix mat(nrows(), ncols() + other.ncols());
		for (size_type i = 0; i < mat.nrows(); i++)
		{
			size_type j;
			for (j = 0; j < ncols(); j++)
			{
				mat(i, j) = data_[i * ncols() + j];
			}
			for (size_type k = 0; k < other.ncols(); k++, j++)
			{
				mat(i, j) = other.data_[i * other.ncols() + k];
			}
		}
		return mat;
	}

	// Iterators and access


	Matrix::reference Matrix::at(size_type i)
	{
		if (i > length())
			throw LinearAlgebraException("Index exceeds flattened matrix length");
		return data_.at(i);
	}

	Matrix::const_reference Matrix::at(size_type i) const
	{
		if (i > length())
			throw LinearAlgebraException("Index exceeds flattened matrix length");
		return data_.at(i);
	}

	Matrix::reference Matrix::at(size_type i, size_type j)
	{
		if (i >= nrows())
			throw LinearAlgebraException("Row index exceeds matrix's size.");
		else if (j >= ncols())
			throw LinearAlgebraException("Column index exceeds matrix's size.");
		return data_.at(j + i * ncols());
	}

	Matrix::const_reference Matrix::at(size_type i, size_type j) const
	{
		if (i >= nrows())
			throw LinearAlgebraException("Row index exceeds matrix's size.");
		else if (j >= ncols())
			throw LinearAlgebraException("Column index exceeds matrix's size.");
		return data_.at(j + i * ncols());
	}

	Matrix::Row Matrix::operator[](size_type i) const
	{
		return Row{ *this, i };
	}

	Matrix::iterator Matrix::begin() { return iterator{ data_.data() }; }
	Matrix::iterator Matrix::end() { return iterator{ data_.data() + data_.size() }; }
	Matrix::const_iterator Matrix::begin() const { return const_iterator{ data_.data() }; }
	Matrix::const_iterator Matrix::end() const { return const_iterator{ data_.data() + data_.size() }; }
	Matrix::row_vector_iterator Matrix::beginRow() { return row_vector_iterator{ *this, 0 }; }
	Matrix::row_vector_iterator Matrix::endRow() { return row_vector_iterator{ *this, nrows() }; }
	Matrix::const_row_vector_iterator Matrix::beginRow() const { return const_row_vector_iterator{ *this, 0 }; }
	Matrix::const_row_vector_iterator Matrix::endRow() const { return const_row_vector_iterator{ *this, nrows() }; }
	Matrix::column_vector_iterator Matrix::beginCol() { return column_vector_iterator{ *this, 0 }; }
	Matrix::column_vector_iterator Matrix::endCol() { return column_vector_iterator{ *this, ncols() }; }
	Matrix::const_column_vector_iterator Matrix::beginCol() const { return const_column_vector_iterator{ *this, 0 }; }
	Matrix::const_column_vector_iterator Matrix::endCol() const { return const_column_vector_iterator{ *this, ncols() }; }

	Matrix::Row Matrix::row(size_type i) { return Row{ *this, i }; }
	Matrix::Column Matrix::col(size_type j) { return Column{ *this, j }; }
	const Matrix::Row Matrix::row(size_type i) const { return Row{ *this, i }; }
	const Matrix::Column Matrix::col(size_type j) const { return Column{ *this, j }; }

	const std::vector<Matrix::value_type>& Matrix::getInternalStdVector() const { return data_; }

	std::string Matrix::to_string() const
	{
		std::ostringstream ostr;
		ostr << *this;
		return ostr.str();
	}

	// Algebraic methods

	Matrix Matrix::inv() const
	{
		size_type m = nrows();
		size_type this_n = ncols();
		if (m != this_n) throw LinearAlgebraException("Matrix is not square.");
		Matrix M = Matrix(*this);
		M = M.concatenate(I(m));
		// gaussian elimination
		size_type n = M.ncols();
		for (size_type i = 0; i < m; i++)
		{ // pivot selection
			size_type pivot = i;
			for (size_type k = i + 1; k < m; k++)
				if (fabs(M(k, i)) > fabs(M(i, i)))
					pivot = k;
			if (pivot > i)
			{ // swapping rows
				for (size_type j = 0; j < n; j++)
				{
					value_type aux = M(i, j);
					M(i, j) = M(pivot, j);
					M(pivot, j) = aux;
				}
			}
			for (size_type k = i + 1; k < m; k++)
			{ // eliminate elements below the pivot
				value_type factor = M(k, i) / M(i, i);
				M(k, i) = static_cast<value_type>(0);
				for (size_type j = i + 1; j < n; j++)
				{
					M(k, j) -= factor * M(i, j);
				}
			}
		}
		// back substitution
		for (size_type j = m - 1; j >= 0 && j < m; j--)
		{ // scale the pivot row to make the pivot element 1
			value_type mjj = M(j, j);
			M(j, j) = 1;
			for (size_type k = j + 1; k < n; k++)
				M(j, k) /= mjj;
			for (size_type i = j - 1; i >= 0 && i < m; i--)
			{ // eliminate elements above the pivot
				value_type mij = M(i, j);
				for (size_type k = j; k < n; k++)
				{
					M(i, k) -= mij * M(j, k);
				}
			}
		}
		M = M.slice(0, m, this_n, this_n * 2);
		return M;
	}


	Matrix::value_type Matrix::determinant_recursion(Matrix M)
	{
		size_type m = M.nrows(), n = M.ncols();
		if (m != n)
			throw LinearAlgebraException("Matrix is not square.");
		if (n == 1)
			return M[0][0];
		if (n == 2)
			return  M[0][0] * M[1][1] - M[0][1] * M[1][0];
		else
		{
			value_type result = 0;
			for (size_type i = 0; i < n; i++)
			{
				Matrix left_submatrix = M.slice(1, n, 0, i);
				Matrix right_submatrix = M.slice(1, n, i + 1, n);
				Matrix submatrix = left_submatrix.concatenate(right_submatrix);
				result += std::pow(-1, i) * M[0][i] * determinant_recursion(submatrix);
			}
			return result;
		}
	}


	Matrix::value_type Matrix::det() const { return determinant_recursion(*this); }


	Matrix Matrix::t() const
	{
		if (isVector())
		{
			if (isRow())
				return Matrix(data_, { ncols(), 1 }); 
			else
				if (isColumn())
					return Matrix(data_, { 1, nrows() });
		}
		else
		{
			size_t m = nrows();
			size_t n = ncols();
			Matrix trans{ n, m };
			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < m; j++)
					trans(i, j) = data_[j * ncols() + i];
			return trans;
		}
		throw LinearAlgebraException("Transpose error: matrix has zero dimensions.");
	}

	Matrix::value_type Matrix::norm(value_type power = static_cast<value_type>(2)) const
	{
		return sqrt(std::accumulate(begin(), end(), static_cast<value_type>(0), [&](value_type accum, value_type next)
			{
				return accum + pow(next, power);
			}));
	}

	Matrix inv(Matrix M) { return M.inv(); }
	Matrix::value_type det(Matrix M) { return M.det(); }
	Matrix t(Matrix M) { return M.t(); }
	Matrix::value_type norm(Matrix M, Matrix::value_type pow) { return M.norm(pow); }

	// Operators

	// Unary operator


	Matrix Matrix::operator-() const
	{
		Matrix mat(getShape());
		std::transform(begin(), end(), mat.begin(), std::negate<value_type>());
		return mat;
	}

	// Operations with scalars


	Matrix Matrix::operator+ (value_type t) const
	{
		Matrix sum{ *this };
		for (reference e : sum) e += t;
		return sum;
	}


	void Matrix::operator+= (value_type t) { for (reference e : data_) e += t; }


	Matrix Matrix::operator- (value_type t) const
	{
		Matrix sum{ *this };
		for (reference e : sum) e -= t;
		return sum;
	}


	void Matrix::operator-= (value_type t) { for (reference e : data_) e -= t; }


	Matrix Matrix::operator* (value_type t) const
	{
		Matrix mul{ *this };
		for (reference e : mul) e *= t;
		return mul;
	}


	void Matrix::operator*= (value_type t) { for (reference e : data_) e *= t; }


	Matrix Matrix::operator/ (value_type t) const
	{
		Matrix mul{ *this };
		for (reference e : mul) e /= t;
		return mul;
	}


	void Matrix::operator/= (value_type t) { for (reference e : data_) e /= t; }

  Matrix operator+ (Matrix::value_type t, const Matrix& M) { return M + t; }
  Matrix operator- (Matrix::value_type t, const Matrix& M) { return M - t; }
  Matrix operator* (Matrix::value_type t, const Matrix& M) { return M * t; }
  Matrix operator/ (Matrix::value_type t, const Matrix& M) { return M / t; }

	// Operations with Vector


	Vector Matrix::operator+(const Vector& v) { return v + *this; }

	void Matrix::operator+=(const Vector& v)
	{
		std::transform(begin(), end(), v.begin(), begin(), std::plus<value_type>());
	}

	Vector Matrix::operator-(Vector v)
	{
		if (!isVector()) throw LinearAlgebraException("Can't subtract Vector from Matrix: Matrix isn't row or column.");
		if (v.size() != length()) throw LinearAlgebraException("Can't subract Vector to Matrix of different length.");
		std::transform(begin(), end(), v.begin(), v.begin(), std::minus<value_type>());
		return v;
	}

	void Matrix::operator-=(const Vector& v)
	{
		std::transform(begin(), end(), v.begin(), begin(), std::minus<value_type>());
	}

	Vector Matrix::operator*(const Vector& v)
	{
		if (v.size() != ncols()) throw LinearAlgebraException("Can't perform dot product: Matrix and Vector have different number of lines.");
		Vector u;
		std::for_each(beginRow(), endRow(), [&](Row row) // make const iterator to allow & argument
			{
				u.push_back(std::inner_product(row.begin(), row.end(), v.begin(), 0.));
			});
		return u;
	}


	// Operations with Matrix


	Matrix Matrix::operator+(const Matrix& other) const
	{
		if (nrows() != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
		if (ncols() != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
		Matrix sum{ getShape() };
		std::transform(begin(), end(), other.begin(), sum.begin(), std::plus<value_type>());
		return sum;
	}


	void Matrix::operator+=(const Matrix& other)
	{
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}


	Matrix Matrix::operator-(const Matrix& other) const
	{
		if (nrows() != other.nrows()) throw LinearAlgebraException("Sum error: matrices have different number of rows.");
		if (ncols() != other.ncols()) throw LinearAlgebraException("Sum error: matrices have different number of columns.");
		Matrix sum{ getShape() };
		std::transform(begin(), end(), other.begin(), sum.begin(), std::minus<value_type>());
		return sum;
	}


	void Matrix::operator-=(const Matrix& other)
	{
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}


	Matrix Matrix::operator*(const Matrix& other) const
	{
		if (other.isScalar())
			return (*this) * other.data_[0];
		else if (isScalar())
			return other * data_[0];
		else if (ncols() != other.nrows()) throw LinearAlgebraException("Multiplication error: matrices do not have appropriate dimensions.");
		else
		{
			Matrix mul(nrows(), other.ncols(), 0);
			std::for_each(beginRow(), endRow(), [&](const Row& row)  // make const iterator to alow & argument
				{
					std::for_each(other.beginCol(), other.endCol(), [&](const Column& col)
						{
							mul.data_.push_back(
								std::inner_product(row.begin(), row.end(), col.begin(), 0.)
							);
						});
				});
			return mul;
		}
	}


	void Matrix::operator*=(const Matrix& other)
	{
		size_type m = nrows();
		if (m != ncols())
			throw LinearAlgebraException("Multiplication error: assigning result of multiplication between non-square matrices to self matrix.");
		if (m != other.nrows() || m != other.ncols())
			throw LinearAlgebraException("Multiplication error: assining result of non-square matrix multiplication to left-hand-side matrix.");
		(*this) = (*this) * other;
	}

	// Operation with Row, Column


	Vector Matrix::operator+(const Row& v) { return v + *this; }

	void Matrix::operator+=(const Row& v)
	{
		std::transform(begin(), end(), v.begin(), begin(), std::plus<value_type>());
	}

	Vector Matrix::operator-(Row r)
	{
		if (!isVector()) throw LinearAlgebraException("Can't subtract from Vector to Matrix: Matrix isn't row or column.");
		if (r.size() != length()) throw LinearAlgebraException("Can't subtract Vector to Matrix of different length.");
		Vector v;
		std::transform(begin(), end(), r.begin(), std::back_inserter(v.getInternalStdVector()), std::minus<value_type>());
		return v;
	}

	void Matrix::operator-=(const Row& v)
	{
		if (!isVector()) throw LinearAlgebraException("Can't subtract from Vector to Matrix: Matrix isn't row or column.");
		if (v.size() != length()) throw LinearAlgebraException("Can't subtract Vector to Matrix of different length.");
		std::transform(begin(), end(), v.begin(), begin(), std::minus<value_type>());
	}

	Vector Matrix::operator*(const Row& v)
	{
		if (v.size() != ncols()) throw LinearAlgebraException("Can't perform product: Matrix and Vector have different number of lines.");
		Vector u;
		std::for_each(beginRow(), endRow(), [&](Row row)  // make const iterator to alow & argument
			{
				u.push_back(std::inner_product(row.begin(), row.end(), v.begin(), 0.));
			});
		return u;
	}


	Vector Matrix::operator+(const Column& v) { return v + *this; }

	void Matrix::operator+=(const Column& v)
	{
		std::transform(begin(), end(), v.begin(), begin(), std::plus<value_type>());
	}

	Vector Matrix::operator-(Column r)
	{
		if (!isVector()) throw LinearAlgebraException("Can't subtract from Vector to Matrix: Matrix isn't row or column.");
		if (r.size() != length()) throw LinearAlgebraException("Can't subtract Vector to Matrix of different length.");
		Vector v;
		std::transform(begin(), end(), r.begin(), std::back_inserter(v.getInternalStdVector()), std::minus<value_type>());
		return v;
	}

	void Matrix::operator-=(const Column& v)
	{
		if (!isVector()) throw LinearAlgebraException("Can't subtract from Vector to Matrix: Matrix isn't row or column.");
		if (v.size() != length()) throw LinearAlgebraException("Can't subtract Vector to Matrix of different length.");
		std::transform(begin(), end(), v.begin(), begin(), std::minus<value_type>());
	}

	Vector Matrix::operator*(const Column& v)
	{
		if (v.size() != ncols()) throw LinearAlgebraException("Can't perform product: Matrix and Vector have different number of lines.");
		Vector u;
		std::for_each(beginRow(), endRow(), [&](Row row) // make const iterator
			{
				u.push_back(std::inner_product(row.begin(), row.end(), v.begin(), 0.));
			});
		return u;
	}


	// Other operations

	std::ostream& operator<<(std::ostream& ostream, const Matrix& matrix)
	{
		ostream << "\n{";
		for (Matrix::size_type i = 0; i < matrix.nrows(); i++)
		{
			ostream << ((i == 0) ? "{ " : " { ");
			for (Matrix::size_type j = 0; j < matrix.ncols(); j++)
				ostream << matrix.at(j + i * matrix.ncols()) << " ";
			ostream << ((i == matrix.nrows() - 1) ? "}" : "}\n");
		}
		ostream << "}\n";
		return ostream;
	}


	void Matrix::operator=(const Vector& v)
	{
		if (!isVector()) throw LinearAlgebraException("Can't assign to Matrix: Matrix isn't row or column.");
		if (length() != v.size()) throw LinearAlgebraException("Can't assign Vector to Matrix of different length.");
		data_.assign(v.begin(), v.end());
	}

	void Matrix::operator=(const Row& v)
	{
		if (!isRow()) throw LinearAlgebraException("Can't assign to Matrix: Matrix isn't row.");
		if (length() != v.size()) throw LinearAlgebraException("Can't assign Row to Matrix of different length.");
		data_.assign(v.begin(), v.end());
	}

	void Matrix::operator=(const Column& v)
	{
		if (!isColumn()) throw LinearAlgebraException("Can't assign to Matrix: Matrix isn't column.");
		if (length() != v.size()) throw LinearAlgebraException("Can't assign Column to Matrix of different length.");
		data_.assign(v.begin(), v.end());
	}

	// END OF MATRIX


	// ROW

	// Iterators and access

	Matrix::iterator Matrix::Row::begin() { return begin_; }
	Matrix::iterator Matrix::Row::end() { return end_; }
	Matrix::const_iterator Matrix::Row::begin() const { return cbegin_; }
	Matrix::const_iterator Matrix::Row::end() const { return cend_; }

	Matrix::size_type Matrix::Row::size() const { return matrix_.ncols(); }

	Matrix::reference Matrix::Row::operator[](size_type i) 
	{ 
		if (i >= size()) throw LinearAlgebraException("Invalid access: index exceeds length of Row.");
		return *(begin_ + i);
	} 

	Matrix::const_reference Matrix::Row::operator[](size_type i) const
	{
		if (i >= size()) throw LinearAlgebraException("Invalid access: index exceeds length of Row.");
		return *(cbegin_ + i);
	}

	std::string Matrix::Row::to_string() const
	{
		std::ostringstream oss;
		oss << *this;
		return oss.str();
	}

		// Algebraic methods

	Matrix::value_type Matrix::Row::norm(value_type power) const
	{
		return sqrt(std::accumulate(begin(), end(), 0., [=](value_type accum, value_type next)
			{
				return accum + pow(next, power);
			}));
	}
	Matrix::value_type norm(Matrix::Row v, Matrix::value_type val) { return v.norm(val); }

	// Operators

	// Unary operator


	Vector Matrix::Row::operator-() const
	{
		Vector v(size());
		std::transform(begin(), end(), v.begin(), std::negate<value_type>());
		return v;
	}

	// Operations with scalar


	Vector Matrix::Row::operator+(value_type t) const
	{
		Vector v(begin(), end());
		for (reference e : v) e += t;
		return v;
	}
	void Matrix::Row::operator+=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e += t; });
	}
	Vector operator+(Matrix::value_type t, const Matrix::Row& v) { return v + t; }


	Vector Matrix::Row::operator-(value_type t) const
	{
		Vector v(begin(), end());
		for (reference e : v) e -= t;
		return v;
	}
	void Matrix::Row::operator-=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e -= t; });
	}
	Vector operator-(Matrix::value_type t, const Matrix::Row& v)
	{
		Vector u(v.begin(), v.end());
		for (Matrix::reference e : u) e = t - e;
		return u;
	}

	Vector Matrix::Row::operator*(value_type t) const
	{
		Vector v(begin(), end());
		for (reference e : v) e *= t;
		return v;
	}

	void Matrix::Row::operator*=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e *= t; });
	}
	Vector operator*(Matrix::value_type t, const Matrix::Row& v) { return v * t; }

	Vector Matrix::Row::operator/(value_type t) const
	{
		Vector v(begin(), end());
		for (reference e : v) e /= t;
		return v;
	}
	void Matrix::Row::operator/=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e /= t; });
	}
	Vector operator/(Matrix::value_type t, const Matrix::Row& v)
	{
		Vector u(v.begin(), v.end());
		for (Matrix::reference e : u) e = t / e;
		return u;
	}

	// Operations with Vector


	Matrix::value_type Matrix::Row::operator*(const Vector& v) const
	{
		return std::inner_product(begin(), end(), v.begin(), static_cast<value_type>(0));
	}

	Vector Matrix::Row::operator+(Vector v) const
	{
		if (size() != v.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), v.begin(), v.begin(), std::plus<value_type>());
		return v;
	}


	void Matrix::Row::operator+=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}


	Vector Matrix::Row::operator-(Vector v) const
	{
		if (size() != v.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), v.begin(), v.begin(), std::minus<value_type>());
		return v;
	}


	void Matrix::Row::operator-=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	// Operations with Matrix

	Vector Matrix::Row::operator+(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Matrix to Row of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}

	void Matrix::Row::operator+=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Row to Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}

	Vector Matrix::Row::operator-(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Matrix of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}

	void Matrix::Row::operator-=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	Vector Matrix::Row::operator*(Matrix M) const
	{
		if (size() != M.nrows()) throw LinearAlgebraException("Can't perform dot product: Row and Matrix have incompatible shapes.");
		Vector v;
		std::for_each(M.beginCol(), M.endCol(), [&](Column col)  // make const iterator to alow & argument
			{
				v.push_back(std::inner_product(col.begin(), col.end(), begin(), static_cast<value_type>(0)));
			});
		return v;
	}

	// Operations with Row


	Matrix::value_type Matrix::Row::operator*(const Row& v) const
	{
		return std::inner_product(begin(), end(), v.begin(), static_cast<value_type>(0));
	}


	Vector Matrix::Row::operator+(const Row& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}


	void Matrix::Row::operator+=(const Row& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}


	Vector Matrix::Row::operator-(const Row& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}


	void Matrix::Row::operator-=(const Row& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	// Operations with Column

	Matrix::value_type Matrix::Row::operator*(const Column& v) const
	{
		return std::inner_product(begin(), end(), v.begin(), static_cast<value_type>(0));
	}

	// Other operations

	std::ostream& operator<<(std::ostream& ostream, const Matrix::Row& row)
	{
		ostream << "{ ";
		for (auto& elem : row) ostream << elem << " ";
		ostream << "}";
		return ostream;
	}

	void Matrix::Row::operator=(const Matrix& m)
	{
		if (!m.isRow()) throw LinearAlgebraException("Can't assign: Matrix isn't row.");
		if (size() != m.length()) throw LinearAlgebraException("Can't assign Matrix to Row of different length.");
		std::copy(m.begin(), m.end(), begin());
	}

	void Matrix::Row::operator=(const Vector& v)
	{
		if (size() != v.size()) throw LinearAlgebraException("Can't assign to Row of different length.");
		std::copy(v.begin(), v.end(), begin());
	}

	void Matrix::Row::operator=(const Row& v)
	{
		if (size() != v.size()) throw LinearAlgebraException("Cannot assign values from Vector to Matrix's Row: different lenghts.");
		std::copy(v.begin(), v.end(), begin());
	}

	// END OF ROW

	// COLUMN

	// Iterators and access

	Matrix::column_iterator Matrix::Column::begin() { return begin_; }
	Matrix::column_iterator Matrix::Column::end() { return end_; }
	Matrix::const_column_iterator Matrix::Column::begin() const { return cbegin_; }
	Matrix::const_column_iterator Matrix::Column::end() const { return cend_; }

	Matrix::size_type Matrix::Column::size() const { return matrix_.nrows(); }

	Matrix::reference Matrix::Column::operator[](size_type i)
	{
		if (i >= size()) throw LinearAlgebraException("Invalid access: index exceeds length of Row.");
		return *(begin_ + i);
	}

	Matrix::const_reference Matrix::Column::operator[](size_type i) const
	{
		if (i >= size()) throw LinearAlgebraException("Invalid access: index exceeds length of Row.");
		return *(cbegin_ + i);
	}

	std::string Matrix::Column::to_string() const
	{
		std::ostringstream oss;
		oss << *this;
		return oss.str();
	}

	// Algebraic methods

	Matrix::value_type Matrix::Column::norm(value_type power) const
	{
		return sqrt(std::accumulate(begin(), end(), 0., [=](value_type accum, value_type next)
			{
				return accum + pow(next, power);
			}));
	}
	Matrix::value_type norm(Matrix::Column v, Matrix::value_type val) { return v.norm(val); }

	// Operators

	// Unary operator


	Vector Matrix::Column::operator-() const
	{
		Vector v(size());
		std::transform(begin(), end(), v.begin(), std::negate<value_type>());
		return v;
	}

	// Operations with scalar


	Vector Matrix::Column::operator+(value_type t) const
	{
		Vector v{ begin(), end() };
		for (reference e : v) e += t;
		return v;
	}
	void Matrix::Column::operator+=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e += t; });
	}
	Vector operator+(Matrix::value_type t, Matrix::Column v) { return v + t; }


	Vector Matrix::Column::operator-(value_type t) const
	{
		Vector v{ begin(), end() };
		for (reference e : v) e -= t;
		return v;
	}
	void Matrix::Column::operator-=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e -= t; });
	}
	Vector operator-(Matrix::value_type t, Matrix::Column v) { return v - t; }

	Vector Matrix::Column::operator*(value_type t) const
	{
		Vector v{ begin(), end() };
		for (reference e : v) e *= t;
		return v;
	}

	void Matrix::Column::operator*=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e *= t; });
	}
	Vector operator*(Matrix::value_type t, Matrix::Column v) { return v * t; }


	Vector Matrix::Column::operator/(value_type t) const
	{
		Vector v{ begin(), end() };
		for (reference e : v) e /= t;
		return v;
	}
	void Matrix::Column::operator/=(value_type t)
	{
		std::for_each(begin(), end(), [&](reference e) { e /= t; });
	}
	Vector operator/(Matrix::value_type t, Matrix::Column v) { return v / t; }

	// Operations with Vector


	Matrix::value_type Matrix::Column::operator*(Vector v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}


	Vector Matrix::Column::operator+(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}


	void Matrix::Column::operator+=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}


	Vector Matrix::Column::operator-(const Vector& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}


	void Matrix::Column::operator-=(const Vector& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	// Operations with Matrix

	Vector Matrix::Column::operator+(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Matrix to Row of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}

	void Matrix::Column::operator+=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Row to Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}

	Vector Matrix::Column::operator-(const Matrix& other) const
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Matrix of different length.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}

	void Matrix::Column::operator-=(const Matrix& other)
	{
		if (!other.isVector()) throw LinearAlgebraException("Can't add to Row: Matrix isn't row or column.");
		if (size() != other.length()) throw LinearAlgebraException("Can't add Matrix of different length.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	Matrix Matrix::Column::operator*(const Matrix& M) const
	{
		if (size() != M.ncols()) throw LinearAlgebraException("Can't perform dot product: Column and Matrix have incompatible shapes.");
		Matrix m{size(), M.ncols()};
		std::for_each(begin(), end(), [&M, &m](const value_type& ele) {
			std::for_each(m.beginRow(), m.endRow(), [&](Row row)  // make const iterator to alow & argument
				{
					row = M * ele;
				});
			});
		return m;
	}

	// Operations with Column


	Matrix::value_type Matrix::Column::operator*(const Column& v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0.);
	}


	Vector Matrix::Column::operator+(const Column& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Columns of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<value_type>());
		return v;
	}


	void Matrix::Column::operator+=(const Column& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Columns of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<value_type>());
	}


	Vector Matrix::Column::operator-(const Column& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Columns of different lengths.");
		Vector v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<value_type>());
		return v;
	}


	void Matrix::Column::operator-=(const Column& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Columns of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<value_type>());
	}

	// Operations with Row

	Matrix Matrix::Column::operator*(Row r) const
	{
		if (size() != r.size()) throw LinearAlgebraException("Can't perform dot product: Column and Matrix have incompatible shapes.");
		Matrix m{ size(), r.size() };
		std::for_each(begin(), end(), [&r, &m](const value_type& ele) {
			std::for_each(m.beginRow(), m.endRow(), [&](Row row)  // make const iterator to alow & argument
				{
					row = r * ele;
				});
			});
		return m;
	}

	// Other operations

	std::ostream& operator<<(std::ostream& ostream, const Matrix::Column& col)
	{
		ostream << "{ ";
		for (auto& elem : col) ostream << elem << " ";
		ostream << "}";
		return ostream;
	}

	void Matrix::Column::operator=(const Matrix& m)
	{
		if (!m.isVector()) throw LinearAlgebraException("Can't assign  Matrix isn't column.");
		if (size() != m.length()) throw LinearAlgebraException("Can't assign Matrix to Column of different length.");
		std::copy(m.begin(), m.end(), begin());
	}

	void Matrix::Column::operator=(const Vector& v)
	{
		if (size() != v.size()) throw LinearAlgebraException("Can't assign to Column of different length.");
		std::copy(v.begin(), v.end(), begin());
	}


} // namespace::alg