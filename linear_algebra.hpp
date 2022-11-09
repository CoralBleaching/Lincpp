#pragma once

#include <string>
#include <vector>
#include <sstream>

/*
	TODO:
	- make containers compatible with STL (and iterators)
	- make equality comparison
	- make const row/col iterator to allow & argument
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
	class Matrix;

	Vector arange(double end, double start = 0., double step = 1.);

	struct Shape {
		size_t m, n;
		bool operator==(const Shape& other) {
			return m == other.m && n == other.n;
		}
	};

	class Iterator;
	class IteratorColumn;
	class Row;
	class Column;
	class IteratorRowVector;
	class IteratorColumnVector;

	class Vector
	{
	private:

		std::vector<double> data_;

	public:
		Vector(size_t k = 0, double val = 0);
		explicit Vector(const std::vector<double> v);
		Vector(const std::initializer_list<double> l);
		Vector(Iterator first, Iterator last);
		Vector(IteratorColumn first, IteratorColumn last);
		explicit Vector(const Matrix& M);


		// Shape methods

		size_t size() const;
		bool isScalar() const;
		void clear();

		Matrix row() const;
		friend Matrix row(Vector v);
		Matrix col() const;
		friend Matrix col(Vector v);
		Matrix t() const;

		// Iterators and access

		std::vector<double>& getInternalStdVector();
		void push_back(double val);
		void insert(Iterator it, double val);
		void insert(Iterator it, const Vector& v);
		Vector concatenate(const Vector&) const;
		Vector slice(size_t start, size_t finish);

		Iterator begin() const;
		Iterator end() const;

		double& operator[](size_t i);
		double& at(size_t i);

		std::string to_string() const;

		// Algebraic methods

		double norm(double power = 2) const;
		friend double norm(Vector v, double val = 2);

		// Operators

		// Unary operation

		Vector operator-() const;

		// Operations with scalar

		Vector operator+(double t) const;
		void operator+=(double t);
		friend Vector operator+(double t, Vector v);

		Vector operator-(double t) const;
		void operator-=(double t);
		friend Vector operator-(double t, Vector v);

		Vector operator*(double t) const;
		void operator*=(double t);
		friend Vector operator*(double t, Vector v);

		Vector operator/(double t) const;
		void operator/=(double t);
		friend Vector operator/(double t, Vector v);

		// Operations with Vector

		double operator*(Vector v) const;
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

		double operator*(const Row& r) const;
		Vector operator+(const Row& r) const;
		Vector operator-(const Row& r) const;
		void operator+=(const Row& r) const;
		void operator-=(const Row& r) const;

		double operator*(const Column& r) const;
		Vector operator+(const Column& r) const;
		Vector operator-(const Column& r) const;
		void operator+=(const Column& r) const;
		void operator-=(const Column& r) const;

		// Other operations

		friend std::ostream& operator<<(std::ostream& ostream, const Vector& vector);

		void operator= (const Matrix&);
		void operator= (const Row&);
		void operator= (const Column&);

	};

	class Matrix
	{
	private:

		static double determinant_recursion(Matrix);

	protected:

		std::vector<double> data_;
		Shape shape_;

		static Shape shapeOf(const std::initializer_list<std::initializer_list<double>>&);

	public: // PRELOAD ALL ROWS AND COLUMNS FOR SPEED EFFICIENCY?
		Matrix(size_t k = 0, Shape shape = { 0, 0 });
		Matrix(size_t m, size_t n);
		Matrix(Shape shape);
		Matrix(std::vector<double> data, Shape shape);
		Matrix(std::vector<double> data, size_t m = -1, size_t n = 1);
		Matrix(Vector vector, Shape shape);
		Matrix(Vector vector, size_t m = -1, size_t n = 1);
		Matrix(const std::initializer_list<std::initializer_list<double>>& nested_list);
		Matrix(const Matrix&) = default;

		// Shape methods

		size_t size() const;
		Shape getShape() const;
		size_t nrows() const;
		size_t ncols() const;
		void setShape(size_t m, size_t n);
		void setShape(Shape shape);
		Matrix reshape(size_t m, size_t n);
		Matrix reshape(Shape shape);
		bool isScalar() const;
		bool isRow() const;
		bool isColumn() const;
		bool isVector() const;
		void clear();

		void push_back(const Vector& v);
		void insert(Iterator it, double val);
		//void insertRows(IteratorRowVector& itIn_beg, IteratorRowVector& itOut_beg, IteratorRowVector& itOut_end);
		//void insertColumns(IteratorColumnVector& itIn_beg, IteratorColumnVector& itOut_beg, IteratorColumnVector& itOut_end);
		Matrix concatenate(const Matrix&) const;
		//Matrix concatenateBelow(const Matrix&);
		Matrix slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end) const;
		//virtual Matrix slice(const arrayIterator& it_b, const arrayIterator& it_e) const;

		// Iterators and access

		Iterator begin() const;
		Iterator end() const;
		IteratorRowVector beginRow() const;
		IteratorRowVector endRow() const;
		IteratorColumnVector beginCol() const;
		IteratorColumnVector endCol() const;

		const double& at(size_t i) const;
		//Row operator[](size_t i) const;
		Row operator[](size_t i) const;
		double& operator()(size_t i, size_t j) { return data_[i * ncols() + j]; }
		Row row(size_t i) const;
		Column col(size_t j) const;

		std::vector<double> getInternalStdVector() const;

		std::string to_string() const;

		/**/
	public:

		// Algebraic methods

		Matrix inv() const;
		friend Matrix inv(Matrix M) { return M.inv(); }
		double det() const;
		friend double det(Matrix M) { return M.det(); }
		Matrix t() const;
		friend Matrix t(Matrix M) { return M.t(); }
		double norm(double) const;
		friend double norm(Matrix M, double pow = 2) { return M.norm(pow); }

		// Operators

		// Unary operator

		Matrix operator- () const;

		// Operations with scalars:

		Matrix operator+ (double) const;
		Matrix operator- (double) const;
		Matrix operator* (double) const;
		Matrix operator/ (double) const;
		void operator+= (double);
		void operator-= (double);
		void operator*= (double);
		void operator/= (double);

		friend Matrix operator+ (double t, const Matrix& M) { return M + t; }
		friend Matrix operator- (double t, const Matrix& M) { return M - t; }
		friend Matrix operator* (double t, const Matrix& M) { return M * t; }
		friend Matrix operator/ (double t, const Matrix& M) { return M / t; }

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

		friend std::ostream& operator<<(std::ostream& ostream, const Matrix& matrix);

		void operator= (const Vector&);
		void operator= (const Row&);
		void operator= (const Column&);
	};

	// MISCELLANEOUS FUNCTIONS

	Matrix t(Vector v);
	Matrix I(size_t);
	Matrix t(Vector v);


	// ITERATORS

	class Iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = double;
		using pointer = double*;
		using reference = double&;
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

		double* getPtr() const { return mptr; }

	private:
		double* mptr;
	};


	class IteratorColumn
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = double;
		using pointer = double*;
		using reference = double&;

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
			return std::distance(it2.getPtr(), it1.getPtr());
		}

		friend bool operator==(const IteratorColumn& a, const IteratorColumn& b) { return a.mptr == b.mptr; }
		friend bool operator!=(const IteratorColumn& a, const IteratorColumn& b) { return a.mptr != b.mptr; }

		double* getPtr() const { return mptr; }
		size_t ncols() const { return ncols_; }

	private:
		double* mptr;
		size_t ncols_;
	};

	class Row
	{
	public:
		Row(Iterator begin = Iterator{ nullptr },
			Iterator end = Iterator{ nullptr },
			Shape shape = { 0,0 }) :
			begin_{ begin }, end_{ end }, shape_{ shape } {} // automatically generate end_?

		Iterator begin() const;
		Iterator end() const;

		Shape getShape() const;
		size_t ncols() const;
		size_t size() const;

		double& operator[](size_t i) const;
		std::string to_string() const;

		// Algebraic methods

		double norm(double power = 2) const;
		friend double norm(Row v, double val = 2);

		// Operators

		// Unary operation

		Vector operator-() const;

		// Operations with scalar

		Vector operator+(double t) const;
		void operator+=(double t);
		friend Vector operator+(double t, Row v);

		Vector operator-(double t) const;
		void operator-=(double t);
		friend Vector operator-(double t, Row v);

		Vector operator*(double t) const;
		void operator*=(double t);
		friend Vector operator*(double t, Row v);

		Vector operator/(double t) const;
		void operator/=(double t);
		friend Vector operator/(double t, Row v);

		// Operations with Vector

		double operator*(const Vector& v) const;
		Vector operator+ (Vector v) const;
		void operator+= (const Vector&);
		Vector operator- (Vector v) const;
		void operator-= (const Vector&);

		// Operations with Matrix

		Vector operator+(const Matrix& M) const;
		Vector operator- (const Matrix&) const;
		void operator+= (const Matrix& M);
		void operator-= (const Matrix& M);
		Vector operator*(Matrix M) const;

		// Operations with Row, Column

		double operator*(const Row&) const;
		Vector operator+ (const Row&) const;
		Vector operator- (const Row&) const;
		void operator+= (const Row&);
		void operator-= (const Row&);
			
		double operator*(const Column&) const;

		// Other operations

		friend std::ostream& operator<<(std::ostream& ostream, const Row& row);

		void operator= (const Matrix&);
		void operator= (const Vector&);
		void operator= (const Row&);

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

		IteratorColumn begin() const;
		IteratorColumn end() const;

		Shape getShape() const;
		size_t nrows() const;
		size_t size() const;

		double& operator[](size_t i);
		std::string to_string() const;

		// Algebraic methods

		double norm(double power = 2) const;
		friend double norm(Column v, double val = 2);

		// Operators

		// Unary operation

		Vector operator-() const;

		// Operations with scalar

		Vector operator+(double t) const;
		void operator+=(double t);
		friend Vector operator+(double t, Column v);

		Vector operator-(double t) const;
		void operator-=(double t);
		friend Vector operator-(double t, Column v);

		Vector operator*(double t) const;
		void operator*=(double t);
		friend Vector operator*(double t, Column v);

		Vector operator/(double t) const;
		void operator/=(double t);
		friend Vector operator/(double t, Column v);

		// Operations with Vector

		double operator*(Vector v) const;
		Vector operator+ (const Vector&) const;
		Vector operator- (const Vector&) const;
		void operator+= (const Vector&);
		void operator-= (const Vector&);

		// Operations with Matrix

		Vector operator+(const Matrix& M) const;
		Vector operator- (const Matrix&) const;
		void operator+= (const Matrix& M);
		void operator-= (const Matrix& M);
		Matrix operator*(const Matrix& M) const;

		// Operations with Row, Column

		double operator*(const Column& v) const;
		Vector operator+ (const Column&) const;
		Vector operator- (const Column&) const;
		void operator+= (const Column&);
		void operator-= (const Column&);

		Matrix operator*(Row r) const;

		// Other operations

		friend std::ostream& operator<<(std::ostream& ostream, const Column& col);

		void operator= (const Matrix&);
		void operator= (const Vector&);
		void operator= (const Column&);

	private:
		IteratorColumn begin_, end_;
		Shape shape_;
	};


	class IteratorRowVector
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = Row;
		using pointer = Iterator;
		using reference = Row&;

		IteratorRowVector(Iterator it = Iterator{ nullptr }, Shape shape = { 0,0 }) : it_{ it }, shape_{ shape } {}
		IteratorRowVector(const IteratorRowVector& rawIterator) = default;
		IteratorRowVector& operator=(const IteratorRowVector& rawIterator) = default;

		Row operator*() { return Row{ it_, it_ + ncols(), shape_ }; }

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

		IteratorColumnVector(const Matrix* ptr = nullptr, size_t index = 0, Shape shape = { 0,0 }) : mptr{ ptr }, index_{ index }, shape_{ shape } {}

		Column operator*() { return mptr->col(index_); }

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

		const Matrix* getMatrixPtr() { return mptr; }
		size_t getIndex() { return index_; }
		size_t nrows() { return shape_.m; }
		size_t ncols() { return shape_.n; }

	private:
		const Matrix* mptr;
		size_t index_;
		Shape shape_;
	};
} // namespace::alg