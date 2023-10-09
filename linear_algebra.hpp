#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <memory>

/*
	TODO:
	- make equality comparison
	- safeguard valid Matrix shapes
	- Vector multiplication shouldn't commute with Column
	- Put bounds on increment/decrement on Iterators (throw)
	- verify const correctness of const_iterators movement methods

	- make containers compatible with STL (and iterators) // DONE
	- make const row/col iterator to allow & argument // DONE
	- const Matrix::row() and const Matrix::col() // DONE
*/

namespace alg {


	class LinearAlgebraException : public std::exception
	{
	public:
		const char* message_;
		LinearAlgebraException(const char* message = "") { message_ = message; }
		virtual const char* what() const throw() { return message_; }
	};

	class Vector;

	struct Shape {
		size_t m, n;
		bool operator==(const Shape& other) {
			return m == other.m && n == other.n;
		}
	};

	class Matrix
	{
	public:
		class Row;
		class Column;
		template<bool> class Iterator;
		template<bool> class ColumnIterator;
		template<bool> class RowVectorIterator;
		template<bool> class ColumnVectorIterator;

	public:
		using value_type = double;
		using size_type = std::size_t;
		using iterator = Iterator<false>;
		using const_iterator = Iterator<true>;
		using column_iterator = ColumnIterator<false>;
		using const_column_iterator = ColumnIterator<true>;
		using row_vector_iterator = RowVectorIterator<false>;
		using const_row_vector_iterator = RowVectorIterator<true>;
		using column_vector_iterator = ColumnVectorIterator<false>;
		using const_column_vector_iterator = ColumnVectorIterator<true>;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type*;
		using reference = value_type&;
		using const_reference = const value_type&;

	private:

		static value_type determinant_recursion(Matrix);

	protected:

		std::vector<value_type> data_;
		Shape shape_;

		static Shape shapeOf(const std::initializer_list<std::initializer_list<value_type>>&);

	public:
		Matrix(Shape shape = { 0, 0 }, value_type val = 0);
		Matrix(size_type m, size_type n, value_type val = 0.);
		Matrix(std::vector<value_type> data, Shape shape);
		Matrix(std::vector<value_type> data, size_type m, size_type n);
		Matrix(std::vector<value_type> data);
		Matrix(Vector vector, Shape shape);
		Matrix(Vector vector, size_type m, size_type n);
		Matrix(Vector vector);
		Matrix(const std::initializer_list<std::initializer_list<value_type>>& nested_list);
		Matrix(const Matrix&) = default;


		// Shape methods

		size_type size() const;
		Shape getShape() const;
		size_type nrows() const;
		size_type ncols() const;
		void setShape(size_type m, size_type n);
		void setShape(Shape shape);
		Matrix reshape(size_type m, size_type n);
		Matrix reshape(Shape shape);
		bool isScalar() const;
		bool isRow() const;
		bool isColumn() const;
		bool isVector() const;
		void clear();

		// Insertions and slice

		void inline insertRow(const Vector& v);
		void inline insertRow(const Vector& v, size_type rowIndex);
		void inline insertRow(const Matrix::Row& v);
		void inline insertRow(const Matrix::Row& v, size_type rowIndex);
		void inline insertRow(const Matrix::Column& v);
		void inline insertRow(const Matrix::Column& v, size_type rowIndex);
		void inline insertRow(const Matrix& v);
		void inline insertRow(const Matrix& v, size_type rowIndex);
		void inline insertRow(const_iterator begin, const_iterator end);
		void inline insertRow(const_column_iterator begin, const_column_iterator end);
		void insertRow(const_iterator begin, const_iterator end, size_type rowIndex);
		void insertRow(const_column_iterator begin, const_column_iterator end, size_type columnIndex);
		void inline insertColumn(const Vector& v);
		void inline insertColumn(const Vector& v, size_type columnIndex);
		void inline insertColumn(const Matrix::Row& v);
		void inline insertColumn(const Matrix::Row& v, size_type columnIndex);
		void inline insertColumn(const Matrix::Column& v);
		void inline insertColumn(const Matrix::Column& v, size_type columnIndex);
		void inline insertColumn(const Matrix& v);
		void inline insertColumn(const Matrix& v, size_type columnIndex);
		void inline insertColumn(const_column_iterator begin, const_column_iterator end);
		void inline insertColumn(const_iterator begin, const_iterator end);
		void insertColumn(const_column_iterator begin, const_column_iterator end, size_type columnIndex);
		void insertColumn(const_iterator begin, const_iterator end, size_type columnIndex);
		Matrix concatenate(const Matrix&) const;
		Matrix slice(size_type row_begin, size_type row_end, size_type column_begin, size_type column_end) const;

		// Iterators and access

		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;
		row_vector_iterator beginRow();
		row_vector_iterator endRow();
		const_row_vector_iterator beginRow() const;
		const_row_vector_iterator endRow() const;
		column_vector_iterator beginCol();
		column_vector_iterator endCol();
		const_column_vector_iterator beginCol() const;
		const_column_vector_iterator endCol() const;

		reference at(size_type i);
		reference at(size_type i, size_type j);
		const_reference at(size_type i) const;
		const_reference at(size_type i, size_type j) const;
		Row operator[](size_type i) const;
		reference operator()(size_type i, size_type j) { return data_[i * ncols() + j]; }
		Row row(size_type i);
		const Row row(size_type i) const;
		Column col(size_type j);
		const Column col(size_type j) const;

		std::vector<value_type> getInternalStdVector() const;

		std::string to_string() const;

		/**/
	public:

		// Algebraic methods

		Matrix inv() const;
		value_type det() const;
		Matrix t() const;
		value_type norm(value_type) const;

		// Operators

		// Unary operator

		Matrix operator- () const;

		// Operations with scalars:

		Matrix operator+ (value_type) const;
		Matrix operator- (value_type) const;
		Matrix operator* (value_type) const;
		Matrix operator/ (value_type) const;
		void operator+= (value_type);
		void operator-= (value_type);
		void operator*= (value_type);
		void operator/= (value_type);

	public:

		// Operations with Vector

		Vector operator+(const Vector& v);
		Vector operator-(Vector v);
		void operator+=(const Vector& v);
		void operator-=(const Vector& v);
		Vector operator*(const Vector& v);

		// Operations with Matrix

		Matrix operator+ (const Matrix&) const;
		void operator+= (const Matrix&);
		Matrix operator- (const Matrix&) const;
		void operator-= (const Matrix&);
		Matrix operator* (const Matrix&) const;
		void operator*= (const Matrix&);

		// Operations with Row, Column

		Vector operator+(const Row& v);
		Vector operator-(Row v);
		void operator+=(const Row& v);
		void operator-=(const Row& v);
		Vector operator*(const Row& v);

		Vector operator+(const Column& v);
		Vector operator-(Column v);
		void operator+=(const Column& v);
		void operator-=(const Column& v);
		Vector operator*(const Column& v);

		// Other operations

		void operator= (const Vector&);
		void operator= (const Row&);
		void operator= (const Column&);
	};

	Matrix inv(Matrix M);
	Matrix::value_type det(Matrix M);
	Matrix t(Matrix M);
	Matrix::value_type norm(Matrix M, Matrix::value_type pow = static_cast<Matrix::value_type>(2));

	Matrix operator+ (Matrix::value_type t, const Matrix& M);
	Matrix operator- (Matrix::value_type t, const Matrix& M);
	Matrix operator* (Matrix::value_type t, const Matrix& M);
	Matrix operator/ (Matrix::value_type t, const Matrix& M);

	std::ostream& operator<<(std::ostream& ostream, const Matrix& matrix);

	// MISCELLANEOUS FUNCTIONS

	Matrix I(size_t);
	Matrix t(Vector v);


	class Vector
	{
	public:
		using value_type = Matrix::value_type;
		using size_type = Matrix::size_type;
		using iterator = Matrix::iterator;
		using const_iterator = Matrix::const_iterator;
		using difference_type = std::ptrdiff_t;
		using reference = Matrix::reference;
		using const_reference = const value_type&;

	private:

		std::vector<value_type> data_;

	public:
		explicit Vector(size_type k = 0, value_type val = 0);
		explicit Vector(const std::vector<value_type> v);
		Vector(const std::initializer_list<value_type> l);
		Vector(const_iterator first, const_iterator last);
		Vector(Matrix::const_column_iterator first, Matrix::const_column_iterator last);
		explicit Vector(const Matrix& M);


		// Shape methods

		size_t size() const;
		bool isScalar() const;
		void clear();

		Matrix row() const;
		Matrix col() const;
		Matrix t() const;

		// Iterators and access

		std::vector<value_type>& getInternalStdVector();
		void push_back(value_type val);
		void insert(iterator it, value_type val);
		void insert(iterator it, const Vector& v);
		Vector concatenate(const Vector&) const;
		Vector slice(size_type start, size_type finish);

		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;

		reference operator[](size_type i);
		reference at(size_type i);

		std::string to_string() const;

		// Algebraic methods

		value_type norm(value_type power = 2) const;

		// Operators

		// Unary operation

		Vector operator-() const;

		// Operations with scalar

		Vector operator+(value_type t) const;
		void operator+=(value_type t);

		Vector operator-(value_type t) const;
		void operator-=(value_type t);

		Vector operator*(value_type t) const;
		void operator*=(value_type t);

		Vector operator/(value_type t) const;
		void operator/=(value_type t);

		// Operations with Vector

		value_type operator*(Vector v) const;
		Vector operator+ (const Vector&) const;
		void operator+= (const Vector&);
		Vector operator- (const Vector&) const;
		void operator-= (const Vector&);

		// Operations with Matrix

		Vector operator+ (const Matrix&) const;
		void operator+= (const Matrix& M);
		Vector operator- (const Matrix&) const;
		void operator-= (const Matrix& M);
		Vector operator*(Matrix M) const;

		// Operations with Row, Column

		value_type operator*(const Matrix::Row& r) const;
		Vector operator+(const Matrix::Row& r) const;
		Vector operator-(const Matrix::Row& r) const;
		void operator+=(const Matrix::Row& r);
		void operator-=(const Matrix::Row& r);

		value_type operator*(const Matrix::Column& r) const;
		Vector operator+(const Matrix::Column& r) const;
		Vector operator-(const Matrix::Column& r) const;
		void operator+=(const Matrix::Column& r);
		void operator-=(const Matrix::Column& r);

		// Other operations


		void operator= (const Matrix&);
		void operator= (const Matrix::Row&);
		void operator= (const Matrix::Column&);

	};

	Matrix row(Vector v);
	Matrix col(Vector v);
	Vector::value_type norm(Vector v, Vector::value_type val = 2);
	Vector operator+(Vector::value_type t, Vector v);
	Vector operator-(Vector::value_type t, Vector v);
	Vector operator*(Vector::value_type t, Vector v);
	Vector operator/(Vector::value_type t, Vector v);
	std::ostream& operator<<(std::ostream& ostream, const Vector& vector);


	Vector arange(
		Matrix::value_type end,
		Matrix::value_type start = static_cast<Matrix::value_type>(0),
		Matrix::value_type step = static_cast<Matrix::value_type>(1)
	);


	// ITERATORS
	
	// Note: the iterator is not trivially constructible.
	template<bool IsConst>
	class Matrix::Iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = difference_type;
		using value_type = value_type;
		using pointer = typename std::conditional_t<IsConst, const value_type*, value_type*>;
		using reference = typename std::conditional_t<IsConst, const value_type&, value_type&>;

		Iterator(pointer ptr) : mptr{ ptr } {}

		template<bool IsOtherConst>
		Iterator(const Iterator<IsOtherConst>& other, std::enable_if_t<IsConst || !IsOtherConst>* = nullptr) : 
			mptr(other.getPtr()) {}

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
		Iterator operator+(difference_type movement) const
		{
			return std::next(*this, movement);
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

		pointer getPtr() const { return mptr; }

	private:
		pointer mptr;
	};

	template<bool IsConst>
	class Matrix::ColumnIterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = value_type;
		using pointer = typename std::conditional_t<IsConst, const value_type*, value_type*>;
		using reference = typename std::conditional_t<IsConst, const value_type&, value_type&>;

		ColumnIterator(pointer ptr, size_type ncols) : ncols_{ ncols } { mptr = ptr; }

		template<bool IsOtherConst>
		ColumnIterator(const Iterator<IsOtherConst>& other, size_type ncols, std::enable_if_t<IsConst || !IsOtherConst>* = nullptr) :
			mptr{ other.getPtr() }, ncols_{ ncols } {}

		template<bool IsOtherConst>
		ColumnIterator(const ColumnIterator<IsOtherConst>& other, std::enable_if_t<IsConst || !IsOtherConst>* = nullptr) :
			mptr{ other.getPtr() }, ncols_{ other.ncols() } {}

		reference operator*() { return *mptr; }
		pointer operator->() { return mptr; }
		ColumnIterator& operator++() { mptr += ncols_; return *this; }
		ColumnIterator operator++(int) { ColumnIterator tmp = *this; ++(*this); return tmp; }
		ColumnIterator& operator--() { mptr -= ncols_; return *this; }
		ColumnIterator operator--(int) { ColumnIterator tmp = *this; --(*this); return tmp; }
		ColumnIterator operator+(difference_type movement)
		{
			auto oldptr = mptr; mptr += movement * ncols_; auto tmp{ *this }; mptr = oldptr;
			return tmp;
		}
		ColumnIterator operator+(difference_type movement) const
		{
			return std::next(*this, movement);
		}
		ColumnIterator operator-(difference_type movement)
		{
			auto oldptr = mptr; mptr -= movement * ncols_; auto tmp{ *this }; mptr = oldptr;
			return tmp;
		}
		void operator+=(difference_type movement) { mptr += movement * ncols_; }
		void operator-=(difference_type movement) { mptr -= movement * ncols_; }
		friend difference_type operator-(const ColumnIterator& it1, const ColumnIterator& it2) // relative distance
		{
			return std::distance(it2.getPtr(), it1.getPtr()) / it1.ncols();
		}
		friend difference_type abs_dist(const ColumnIterator& it1, const ColumnIterator& it2) // absolute distance
		{
			return std::distance(it2.getPtr(), it1.getPtr());
		}

		friend bool operator==(const ColumnIterator& a, const ColumnIterator& b) { return a.mptr == b.mptr; }
		friend bool operator!=(const ColumnIterator& a, const ColumnIterator& b) { return a.mptr != b.mptr; }
		pointer getPtr() const { return mptr; }
		size_type ncols() const { return ncols_; }

	private:
		pointer mptr;
		size_type ncols_;
	};

	class Matrix::Row
	{
	public:
		Row(const Matrix& matrix, size_type row_index) :
			matrix_{ matrix },
			row_index_{ row_index },
			begin_{ const_cast<pointer>((matrix.begin() + row_index * matrix.ncols()).getPtr()) },
			end_{ const_cast<pointer>((matrix.begin() + (row_index + 1) * matrix.ncols()).getPtr()) },
			cbegin_{ matrix.begin() + row_index * matrix.ncols() },
			cend_{ matrix.begin() + (row_index + 1) * matrix.ncols() }
		{}

		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;

		size_type size() const;

		reference operator[](size_type i);
		const_reference operator[](size_type i) const;

		std::string to_string() const;

		// Algebraic methods

		value_type norm(value_type power = static_cast<value_type>(2)) const;

		// Operators

		// Unary operation

		Vector operator-() const;

		// Operations with scalar

		Vector operator+(value_type t) const;
		void operator+=(value_type t);
		void operator+=(value_type) const = delete;
		friend Vector operator+(value_type t, const Row& v);

		Vector operator-(value_type t) const;
		void operator-=(value_type t);
		void operator-=(value_type t) const = delete;
		friend Vector operator-(value_type t, const Row& v);

		Vector operator*(value_type t) const;
		void operator*=(value_type t);
		void operator*=(value_type t) const = delete;
		friend Vector operator*(value_type t, const Row& v);

		Vector operator/(value_type t) const;
		void operator/=(value_type t);
		void operator/=(value_type t) const = delete;
		friend Vector operator/(value_type t, const Row& v);

		// Operations with Vector

		value_type operator*(const Vector& v) const;
		Vector operator+ (Vector v) const;
		void operator+= (const Vector&);
		void operator+= (const Vector&) const = delete;
		Vector operator- (Vector v) const;
		void operator-= (const Vector&);
		void operator-= (const Vector&) const = delete;

		// Operations with Matrix

		Vector operator+(const Matrix& M) const;
		Vector operator- (const Matrix&) const;
		void operator+= (const Matrix& M);
		void operator+= (const Matrix& M) const = delete;
		void operator-= (const Matrix& M);
		void operator-= (const Matrix& M) const = delete;
		Vector operator*(Matrix M) const;

		// Operations with Row, Column

		value_type operator*(const Row&) const;
		Vector operator+ (const Row&) const;
		Vector operator- (const Row&) const;
		void operator+= (const Row&);
		void operator+= (const Row&) const = delete;
		void operator-= (const Row&);
		void operator-= (const Row&) const = delete;
			
		value_type operator*(const Column&) const;

		// Other operations

		friend std::ostream& operator<<(std::ostream& ostream, const Row& row);

		void operator= (const Matrix&);
		void operator= (const Vector&);
		void operator= (const Row&);
		void operator= (const Matrix&) const = delete;
		void operator= (const Vector&) const = delete;
		void operator= (const Row&) const = delete;

	protected:
		const Matrix& matrix_;
		size_type row_index_;
		iterator begin_, end_;
		const_iterator cbegin_, cend_;
	};

	Matrix::value_type norm(Matrix::Row v, Matrix::value_type val = static_cast<Matrix::value_type>(2));


	class Matrix::Column
	{
	public:
		Column(const Matrix& matrix, size_type col_index) :
			matrix_{ matrix }, col_index_{ col_index },
			begin_{ const_cast<pointer>((matrix.begin() + col_index).getPtr()), matrix.ncols()},
			end_{ const_cast<pointer>((matrix.begin() + matrix.ncols() * matrix.nrows() + col_index).getPtr()), matrix.ncols()},
			cbegin_{ matrix.begin() + col_index, matrix.ncols() },
			cend_{ matrix.begin() + matrix.ncols() * matrix.nrows() + col_index, matrix.ncols() }
		{}

		//begin() + j,
		//begin() + ncols() * nrows() + j,

		column_iterator begin();
		column_iterator end();
		const_column_iterator begin() const;
		const_column_iterator end() const;

		size_type size() const;

		reference operator[](size_type i);
		const_reference operator[](size_type i) const;
		std::string to_string() const;

		// Algebraic methods

		value_type norm(value_type power = static_cast<value_type>(2)) const;

		// Operators

		// Unary operation

		Vector operator-() const;

		// Operations with scalar

		Vector operator+(value_type t) const;
		void operator+=(value_type t);
		void operator+=(value_type t) const = delete;
		friend Vector operator+(value_type t, Column v);

		Vector operator-(value_type t) const;
		void operator-=(value_type t);
		void operator-=(value_type t) const = delete;
		friend Vector operator-(value_type t, Column v);

		Vector operator*(value_type t) const;
		void operator*=(value_type t);
		void operator*=(value_type t) const = delete;
		friend Vector operator*(value_type t, Column v);

		Vector operator/(value_type t) const;
		void operator/=(value_type t);
		void operator/=(value_type t) const = delete;
		friend Vector operator/(value_type t, Column v);

		// Operations with Vector

		value_type operator*(Vector v) const;
		Vector operator+ (const Vector&) const;
		Vector operator- (const Vector&) const;
		void operator+= (const Vector&);
		void operator+= (const Vector&) const = delete;
		void operator-= (const Vector&);
		void operator-= (const Vector&) const = delete;

		// Operations with Matrix

		Vector operator+(const Matrix& M) const;
		Vector operator-(const Matrix& M) const;
		void operator+= (const Matrix& M);
		void operator+= (const Matrix& M) const = delete;
		void operator-= (const Matrix& M);
		void operator-= (const Matrix& M) const = delete;
		Matrix operator*(const Matrix& M) const;

		// Operations with Row, Column

		value_type operator*(const Column& v) const;
		Vector operator+ (const Column&) const;
		Vector operator- (const Column&) const;
		void operator+= (const Column&);
		void operator+= (const Column&) const = delete;
		void operator-= (const Column&);
		void operator-= (const Column&) const = delete;

		Matrix operator*(Row r) const;

		// Other operations

		friend std::ostream& operator<<(std::ostream& ostream, const Column& col);

		void operator= (const Matrix&);
		void operator= (const Vector&);
		void operator= (const Column&);
		void operator= (const Matrix&) const = delete;
		void operator= (const Vector&) const = delete;
		void operator= (const Column&) const = delete;

	private:
		const Matrix& matrix_;
		size_type col_index_;
		column_iterator begin_, end_;
		const_column_iterator cbegin_, cend_;
	};

	Matrix::value_type norm(Matrix::Column v, Matrix::value_type val = static_cast<Matrix::value_type>(2));

	template<bool IsConst>
	class Matrix::RowVectorIterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = Row;
		using pointer = typename std::conditional_t<IsConst, const Row*, Row*>;
		using reference = typename std::conditional_t<IsConst, const Row&, Row&>;
		using matrix_reference = typename std::conditional_t<IsConst, const Matrix&, Matrix&>;

		RowVectorIterator(matrix_reference matrix, size_type row_index) : matrix_{ matrix }, row_index_{ row_index } {}

		template<bool IsOtherConst>
		RowVectorIterator(const RowVectorIterator<IsOtherConst>& other, std::enable_if_t<IsConst || !IsOtherConst>* = nullptr) :
			matrix_{ other.matrix_ }, row_index_{ other.row_index_ } {}

		value_type operator*() { return value_type(matrix_, row_index_); }
		//pointer operator->() { return std::addressof(matrix_.row(row_index_)); }

		RowVectorIterator& operator++() { ++row_index_; return *this; }
		RowVectorIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		RowVectorIterator& operator--() { --row_index_; return *this; }
		RowVectorIterator operator--(int) { auto tmp = *this; --(*this); return tmp; }
		RowVectorIterator operator+(difference_type movement)
		{
			auto oldidx = row_index_; row_index_ += movement; auto tmp{ *this }; row_index_ = oldidx;
			return tmp;
		}
		RowVectorIterator operator+(difference_type movement) const
		{
			return std::next(*this, movement);
		}
		RowVectorIterator operator-(difference_type movement)
		{
			auto oldidx = row_index_; row_index_ -= movement; auto tmp{ *this }; row_index_ = oldidx;
			return tmp;
		}
		RowVectorIterator operator-(difference_type movement) const
		{
			return std::prev(*this, movement);
		}
		void operator+=(difference_type movement) { row_index_ += movement; }
		void operator-=(difference_type movement) { row_index_ -= movement; }
		friend bool operator==(const RowVectorIterator& a, const RowVectorIterator& b) 
		{ 
			return a.row_index_ == b.row_index_ && std::addressof(a.matrix_) == std::addressof(b.matrix_);
		}
		friend bool operator!=(const RowVectorIterator& a, const RowVectorIterator& b) { return !(a == b); }
		inline size_type getIndex() const { return row_index_; }

	private:
		matrix_reference matrix_;
		size_type row_index_;
	};

	template<bool IsConst>
	class Matrix::ColumnVectorIterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = Column;
		using pointer = typename std::conditional_t<IsConst, const Column*, Column*>;
		using reference = typename std::conditional_t<IsConst, const Column&, Column&>;
		using matrix_reference = typename std::conditional_t<IsConst, const Matrix&, Matrix&>;

		ColumnVectorIterator(matrix_reference matrix, size_type col_index) : matrix_{ matrix }, col_index_{ col_index } {}

		template<bool IsOtherConst>
		ColumnVectorIterator(const ColumnVectorIterator<IsOtherConst>& other, std::enable_if_t<IsConst || !IsOtherConst>* = nullptr) :
			matrix_{ other.matrix_ }, col_index_{ other.col_index_ } {}

		value_type operator*() { return value_type(matrix_, col_index_); }
		//pointer operator->() { return std::addressof(matrix_.col(col_index_)); }

		ColumnVectorIterator& operator++() { col_index_ += 1; return *this; }
		ColumnVectorIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		ColumnVectorIterator& operator--() { col_index_ -= 1; return *this; }
		ColumnVectorIterator operator--(int) { auto tmp = *this; --(*this); return tmp; }
		ColumnVectorIterator operator+(difference_type movement)
		{
			auto oldidx = col_index_; col_index_ += movement; auto tmp{ *this }; col_index_ = oldidx;
			return tmp;
		}
		ColumnVectorIterator operator-(difference_type movement)
		{
			auto oldidx = col_index_; col_index_ -= movement; auto tmp{ *this }; col_index_ = oldidx;
			return tmp;
		}
		void operator+=(difference_type movement) { col_index_ += movement; }
		void operator-=(difference_type movement) { col_index_ -= movement; }

		friend bool operator==(const ColumnVectorIterator& a, const ColumnVectorIterator& b)
		{
			return a.col_index_ == b.col_index_ && std::addressof(a.matrix_) == std::addressof(b.matrix_);
		}
		friend bool operator!=(const ColumnVectorIterator& a, const ColumnVectorIterator& b) { return !(a == b); }
		inline size_type getIndex() const { return col_index_; }

	private:
		matrix_reference matrix_;
		size_type col_index_;
	};
} // namespace::alg