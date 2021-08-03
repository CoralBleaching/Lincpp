#include<vector>
#include<sstream>
#include<ostream>
#include<cmath>


template <typename T>
std::ostream& operator<< (std::ostream&, std::vector<T>);

namespace alg {

	template<class T>
	class Matrix;

	/*
		Exception type to specify that an error/exception happened within the linear algebra routines below.
		The message is customizable to inform about the exact nature of the linear algebra error.
	*/
	class linear_algebra_exception : public std::exception
	{
	public:

		const char* message;

		linear_algebra_exception(const char* message = "")
		{
			this->message = message;
		}

		virtual const char* what()  const throw()
		{
			return this->message;
		}
	};

	/*
		A type that encapsulates the behavior of mathematical vectors as defined in linear algebra treatments. 
		It is a wrapper around std::vector defining mathematical operations on it in a proper namespace (for evading conflicts).
		Has 5 constructors:
			default           vector();
			copy              vector(v);
			matrix            vector(Matrix);
			standard          vector(n); 
			initializer_list  vector({ v1, v2, ..., vn });
		Supports the most principal std::vector's methods: begin, end, push_back, size.
		Supports brackets and parenthesis indexing.
		alg::vector objects can be printed to standard output.
		Mathematical perators on alg::vector are defined to work as mathematical operations.
		Example:
			vector<double> v1 = { 1, 2, 3 };
			vector<double> v2 = 1 + v1;      ----> v2: { 2, 3, 4 }
			double p = v1 * v2;              ----> p: 20 
		Note: alg::vector does not distinguish between row or column vectors! In order to make that distinction,
		it is recommended to use an alg::Matrix object.
	*/
	template<class T>
	class vector 
	{
		
	private:
		
		std::vector<T> _vector; // Underlying data structure is a standard vector.

	public:

		// Constructors

		// Default constructor
		vector()
		{
			_vector = std::vector<T>();
		}

		// Copy constructor
		vector(std::vector<T> user_v)
		{
			_vector = user_v;
		}

		// Matrix constructor
		vector(Matrix<T> user_m)
		{
			try
			{
				*this = user_m.toVector();
			}
			catch (linear_algebra_exception e)
			{
				throw e;
			}
		}

		// Standard allocating constructor. Initializes the entries with a given value, 0 by default.
		vector(size_t m, T value = 0)
		{
			_vector = std::vector<T>(m, value);
		}

		// Constructor with initializer list
		vector(std::initializer_list<T> list) : vector((std::vector<T>)list) {};
		
		struct Iterator
		{
			using iterator_category = std::bidirectional_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using value_type = T;
			using pointer = T*;
			using reference = T&;

			Iterator(pointer ptr) : m_ptr(ptr) {}

			reference operator*() { return *m_ptr; }
			pointer operator->() { return m_ptr; }
			Iterator& operator++() { m_ptr++; return *this; }
			Iterator& operator++(int) { Iterator temp = *this; ++(*this); return temp; }
			Iterator& operator--() { m_ptr--; return *this; }
			Iterator& operator--(int) { Iterator temp = *this; --(*this); return temp; }
			friend bool operator==(const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; }
			friend bool operator!=(const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; }

			pointer m_ptr;
		};

		// Adapted methods from std::vector

		Iterator begin()
		{
			return Iterator(_vector.begin()._Ptr);
		}

		Iterator end()
		{
			return Iterator(_vector.end()._Ptr);
		}
 
		size_t size()
		{
			return _vector.size();
		}

		void push_back(T t)
		{
			_vector.push_back(t);
		}

		void pop_back()
		{
			_vector.pop_back();
		}

		T& back()
		{
			return *(--end());
		}
		/**/
		void erase(Iterator i)
		{
			vector<T> new_vector = {};
			for (auto it = begin(); it != end(); it++)
			{
				if (it != i)
					new_vector.push_back(*it);
			}
			list = new_vector;
		}
 
		void insert(Iterator it1, Iterator it2_b, Iterator it2_e)
		{
			std::vector<double> new_vector = {};
			for (auto it = begin(); it != ++end(); it++)
			{
				if (it == it1)
				{
					for (auto jt = it2_b; jt != it2_e; jt++)
						new_vector.push_back(*jt);
				}
				if (it != end())
					new_vector.push_back(*it);
			}
			_vector = new_vector;
		}
 
		// Element access operators 

		T& operator[] (size_t i)
		{
			return _vector[i];
		}

		T& operator() (size_t i)
		{
			return _vector[i];
		}

		// Assignment operator

		void operator=(vector<T> other)
		{
			_vector = other._vector;
		}

		void operator=(std::initializer_list<T> list)
		{
			_vector = std::vector<T>(list);
		}

		// Comparison operators

		friend bool operator==(const vector<T>& lhs, const vector<T>& rhs)
		{
			return lhs._vector == rhs._vector;
		}

		friend bool operator!=(const vector<T>& lhs, const vector<T>& rhs)
		{
			return lhs._vector != rhs._vector;
		}
		/*/
		// Scalar operator. 
		// If a vector has size 1, it can be (automatically) cast into its scalar type.
		operator T()
		{
			if (!isScalar())
				throw linear_algebra_exception("Tried to convert a non-scalar vector into a T.");
			return _vector[0];
		}
		/**/
		// Returns true if the number of elements in vector is 1, false otherwise.
		bool isScalar()
		{
			return size() == 1;
		}

		// Concatenates the elements of another vector to the right, returning a new, longer vector.
		vector<T> concatenate(vector<T> other)
		{
			vector v = *this;
			size_t m = size();
			v.insert(v.end(), other.begin(), other.end());
			return v;
		}

		// Returns a vector comprised of a sublist of elements within a range of indexes of the original vector.
		vector<T> slice(size_t index_begin, size_t index_end)
		{
			size_t m = index_end - index_begin;
			if (index_end > size())
				throw linear_algebra_exception("Slice exceeded vector dimension.");
			vector v = vector(m);
			for (size_t i = 0; i < m; i++)
				v[i] = _vector[i + index_begin];
			return v;
		}

		// The Euclidean norm of the vector.
		T norm()
		{
			T sum = 0;
			size_t m = size();
			for (size_t i = 0; i < m; i++)
				sum += _vector[i] * _vector[i];
			return std::sqrt(sum);
		}

		// The Euclidean norm of the vector v. (Alternative for clarity of sintax)
		friend T norm(vector<T> v)
		{
			return v.norm();
		}

		// Unary operator

		vector<T> operator- ()
		{
			vector v = *this;
			for (auto& vii : v)
				vii = -vii;
			return v;
		}

		// Operations with scalars:

		vector<T> operator+ (T t)
		{
			vector v = *this;
			size_t m = size();
			for (size_t i = 0; i < m; i++)
				v[i] += t;
			return v;
		}

		friend vector<T> operator+ (T t, vector<T> v)
		{
			return v + t;
		}

		vector<T> operator- (T t)
		{
			return v + (-t);
		}

		friend vector<T> operator- (T t, vector<T> v)
		{
			return v - t;
		}

		vector<T> operator* (T t)
		{
			vector v = *this;
			size_t m = size();
			for (size_t i = 0; i < m; i++)
				v[i] *= t;
			return v;
		}

		friend vector<T> operator* (T t, vector<T> v)
		{
			return v * t;
		}

		vector<T> operator/ (T t)
		{
			vector v = *this;
			size_t m = size();
			for (size_t i = 0; i < m; i++)
				v[i] /= t;
			return v;
		}

		friend vector<T> operator/ (T t, vector<T> v)
		{
			return v / t;
		}

		// Operations with alg::vector

		vector<T> operator+ (vector<T> other)
		{
			size_t m = size();
			if (m != other.size()) throw linear_algebra_exception("Sum error: vectors have different number of columns.");
			vector sum = vector(m);
			for (size_t i = 0; i < m; i++)
				sum[i] = _vector[i] + other[i];
			return sum;
		}

		void operator+= (vector<T> other)
		{
			*this = *this + other;
		}

		vector<T> operator- (vector<T> other)
		{
			return *this + (-other);
		}

		void operator-= (vector<T> other)
		{
			*this = *this - other;
		}

		T operator* (vector<T> other)
		{
			size_t m = size();
			if (m != other.size())
				throw linear_algebra_exception("Multiplication error: vectors do not have equal dimensions.");
			else
			{
				T mul = 0;
				for (size_t i = 0; i < m; i++)
					mul += _vector[i] * other[i];
				return mul;
			}
		}

		// Overload allows for easy printing of vector's contents.
		friend std::ostream& operator<< (std::ostream& output, vector<T> v)
		{
			size_t m = v.size();
			std::ostringstream line;
			line << "{ ";
			for (size_t j = 0; j < m; j++)
			{
				if (j < m - 1) line << v[j] << ", ";
				else line << v[j] << " }";
			}
			output << line.str();
			return output;
		}
	};

	/*
		A type that encapsulates the behavior of mathematical matrices as defined in linear algebra treatments.
		It is built upon the alg::vector class defined above. It is a wrapper around a container of alg::vectors 
		defining mathematical operations on it. 
		Note: there's no safeguard against violation the assumptions about a matrix's shape (it doesn't check whether
			  different lines in a single matrix all have the same length).
		Has 6 constructors:
			default           Matrix();
			copy              Matrix(M);
			matrix            Matrix(Matrix);
			standard          Matrix(n);
			std::vector       Matrix(std::vector<alg::vector>)
			initializer_list  vector({ { a_11, ..., a_1n } , ... , { a_n1, ..., a_nn } });
		Supports the following std::vector methods: begin, end, push_back, size.
		Supports brackets and parenthesis indexing.
		Can be printed to standard output.
		Example:
			Matrix<double> m = { {1, 2}, {3, 4} };
			cout << m(1,0);  // ----> 3
			cout << m[0];    // ----> { 1, 2 }
			cout << m[0][1]; // ----> 2
			cout << m;       // ----> { {1, 2}, {3, 4} };
		Mathematical operators on alg::vector are defined to work as mathematical operations. The philosophy of 
		operations is to always "downgrade" into the lesser dimensional type after operations with lesser dimensional
		types. This should greatly simplify sintax in the code.
		Example:
			Matrix<double> m = { {1, 2}, {3, 4} };
			vector<double> v = { 5, 6 };
			auto p = m * v;   
			cout << p;      //   ----> { 17, 39 } - type: alg::vector
		Methods available include:
			- determinant and norm
			- t (transpose)
			- inverse
			- concatenate and slice
			- column and row
			- toVector
		TBD:
			- eigenvalues and eigenvectors
			- factorization
	*/
	template <class T>
	class Matrix {

	private:

		std::vector<vector<T>> _matrix; // Underlying data structure is a std::vector of alg::vectors.

	public:

		// Constructors

		// Default constructor
		Matrix()
		{
			_matrix = std::vector<vector<T>>();
		}

		// vector constructor. Creates a (matrix object encapsulating a) column vector by default.
		Matrix(vector<T> user_v, bool column = true)
		{
			if (column)
			{
				size_t m = user_v.size();
				_matrix = std::vector<vector<T>>(m);
				for (size_t i = 0; i < m; i++) (*this)[i] = { user_v[i] };
			}
			else
				_matrix = { user_v };
		}

		// Copy constructor
		Matrix(Matrix<T>& user_m)
		{
			size_t m = user_m.n_rows();
			_matrix = std::vector<vector<T>>(m);
			for (size_t i = 0; i < m; i++) _matrix[i] = user_m[i];
		}

		// Standard allocating constructor
		// Initializes all entries (if applicable) with (T)0.
		Matrix(size_t m, size_t n = 0)
		{
			_matrix = std::vector<vector<T>>(m);
			for (size_t i = 0; n != 0 && i < m; i++)
			{
				_matrix[i] = vector<T>(n);
			}
		}

		// std::vector constructors (necessary for initializer_list constructor implementation)
		Matrix(std::vector<vector<T>> user_m)
		{
			size_t m = user_m.size();
			_matrix = std::vector<vector<T>>(m);
			for (size_t i = 0; i < m; i++)
				_matrix[i] = user_m[i];
		}

		// Initializer list constructor
		Matrix(std::initializer_list<std::initializer_list<T>> list)
		{
			*this = Matrix();
			for (auto& sublist : list)
				_matrix.push_back((vector<T>)sublist);
		}
		
		// Methods inherited from std::vector

		auto begin()
		{
			return _matrix.begin();
		}

		auto end()
		{
			return _matrix.end();
		}

		void push_back(vector<T> v)
		{
			if (_matrix.size() > 0 && v.size() != n_columns())
				throw linear_algebra_exception("Error: pushing back vector with inappropriate dimensions.");
			_matrix.push_back(v);
		}

		// Element access operators

		vector<T>& operator[] (size_t i)
		{
			return _matrix[i];
		}

		vector<T>& operator() (size_t i)
		{
			return _matrix[i];
		}

		T& operator() (size_t i, size_t j)
		{
			return _matrix[i][j];
		}

		// Returns the number of rows
		size_t n_rows()
		{
			return _matrix.size();
		}


		// Returns the number of columns
		size_t n_columns()
		{
			return _matrix[0].size();
		}

		// Returns row m of the matrix as an alg::vector. m starts from 0. 
		vector<T> row(size_t m)
		{
			if (n_rows() < m) throw linear_algebra_exception("Access error: trying to access row out of index.");
			return vector<T>(_matrix[m]);
		}

		/* 
			Returns column n of the matrix as (an alg::Matrix representing) a column vector, i.e. a 1xn matrix.
			n starts from 0.
		*/
		Matrix<T> column(size_t n)
		{
			if (n_columns() < n) throw linear_algebra_exception("Access error: trying to access column out of index.");
			size_t m = n_rows();
			Matrix<T> col(m);
			for (size_t i = 0; i < m; i++) col[i].push_back(_matrix[i][n]);
			return col;
		}

		// Concatenates the elements of another Matrix to the right. Both matrices must have the same number of rows.
		Matrix<T> concatenate(Matrix<T> other)
		{
			Matrix M = Matrix(*this);
			size_t m = n_rows();
			if (m != other.n_rows()) throw linear_algebra_exception("Matrices of different number of rows.");
			for (size_t i = 0; i < m; i++)
				M[i].insert(M[i].end(), other[i].begin(), other[i].end());
			return M;
		}

		/*
			Creates a submatrix of this matrix with the elements ranging from row 'row_begin' to 'row_end - 1' 
			and column 'column_begin' to 'column_end - 1', then returns it.
		*/
		Matrix<T> slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end)
		{
			size_t m = row_end - row_begin;
			size_t n = column_end - column_begin;
			if (row_end > n_rows() || column_end > n_columns()) throw linear_algebra_exception("Slice exceeded matrix dimensions.");
			Matrix M = Matrix(m, n);
			for (size_t i = 0; i < m; i++)
				for (size_t j = 0; j < n; j++)
					M[i][j] = _matrix[i + row_begin][j + column_begin];
			return M;
		}

		/*
			Returns a new Matrix object containing the inverse of the current matrix. It applies 
			an algorithm for LU decomposition.
		*/
		Matrix<T> inverse()
		{
			size_t m = n_rows();
			size_t this_n = n_columns();
			if (m != this_n) throw linear_algebra_exception("Matrix is not square.");
			Matrix M = Matrix<T>(*this);
			M = M.concatenate(I<T>(m));
			size_t n = M.n_columns();
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

		/*
			Internal routine from calculating the determinant of matrix M via the cofactor formula.
			INEFFICIENT, but coded in a hurry.
		*/
		T determinant_recursion(Matrix M)
		{
			size_t m = M.n_rows(), n = M.n_columns();
			if (m != n)
				throw linear_algebra_exception("Matrix is not square.");
			if (n == 1)
				return M[0][0];
			if (n == 2)
				return M[0][0] * M[1][1] - M[0][1] * M[1][0];
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

		// Main method returning the determinant of the matrix. Calls a recursive subroutine on itself.
		T determinant()
		{
			return determinant_recursion(*this);
		}

		// Returns true if the matrix is a singleton, i.e. it has only one entry. Returns false otherwise.
		bool isScalar()
		{
			return n_columns() == 1 && n_rows() == 1;
		}

		// Returns true if the matrix has a single column, false otherwise.
		bool isColumn()
		{
			return n_columns() == 1;
		}

		// Returns true if the matrix has a single row, false otherwise. 
		bool isRow()
		{
			return this->n_rows() == 1;
		}

		// Returns true if the matrix has either a single column or a single row, i.e. it's a representation of a mathematical vector.
		bool isVector()
		{
			return this->isRow() || this->isColumn();
		}

		// If the matrix has either a single column or a single row, converts it into an alg::vector (erasing information about orientation).
		vector<T> toVector()
		{
			if (isRow())
				return vector<T>(_matrix[0]);
			else if (isColumn())
			{
				vector<T> v = {};
				size_t m = n_rows();
				for (size_t i = 0; i < m; i++)
					v.push_back(_matrix[i][0]);
				return v;
			}
	//			return vector<T>(column(0));
			else
				throw linear_algebra_exception("Cannot convert non-row/non-column matrix to vector.");
		}

		// Returns a transposed copy of matrix M.
		friend Matrix<T> t(Matrix<T> M)
		{
			if (M.isVector())
			{
				if (M.isRow())
					return Matrix(M[0], true); // this constructor creates a column matrix.
				else
					if (M.isColumn())
					{
						size_t n = M.n_rows();
						Matrix matrix = Matrix<T>(1, n);
						for (size_t j = 0; j < n; j++) matrix[0][j] = M[j][0];
						return matrix;
					}
			}
			else
			{
				Matrix matrix = Matrix<T>(M);
				size_t m = M.n_rows();
				size_t n = M.n_columns();
				if (m != n)
				{
					std::ostringstream error_msg("Transpose error: non-square matrix. Dimensions: ");
					error_msg << m << "x" << n << std::endl;
					throw linear_algebra_exception(error_msg.str().data());
				}
				for (size_t i = 0; i < m; i++)
					for (size_t j = 0; j < i; j++)
					{
						T aux = matrix[i][j];
						matrix[i][j] = matrix[j][i];
						matrix[j][i] = aux;
					}
				return matrix;
			}
			throw linear_algebra_exception("Transpose error: matrix has zero dimensions.");
		}

		// Method for creating a column vector (of alg::Matrix type) from an alg::vector object.
		friend Matrix<T> t(vector<T> v)
		{
			Matrix<T> M = Matrix<T>(v);
			return M;
		}

		// The Euclidean norm of the row/column matrix.
		T norm()
		{
			if (isVector())
			{
				vector<T> v = toVector();
				return v.norm();
			}
			throw linear_algebra_exception("Tried to take the norm of non-row/non-column matrix.");
		}

		// The Euclidean norm of the row/column matrix. (Alternative for clarity of sintax)
		friend T norm(Matrix<T> v)
		{
			return v.norm();
		}

		// Overload allows for easy printing of matrix's contents.
		friend std::ostream& operator<< (std::ostream& output, Matrix<T> matrix)
		{
			size_t m = matrix.n_rows();
			size_t n = matrix.n_columns();
			std::ostringstream line;
			line << "{";
			for (size_t i = 0; i < m; i++)
			{
				std::ostringstream column;
				column << "{ ";
				for (size_t j = 0; j < n; j++)
				{
					if (j < n - 1) column << matrix[i][j] << ", ";
					else column << matrix[i][j] << " }";
				}
				if (i > 0) line << " ";
				if (i < m - 1) line << column.str() << ",\n";
				else line << column.str() << "}";
			}
			output << line.str();
			return output;
		}

		// Mathematical operators overloading

		// Unary operator

		Matrix<T> operator- ()
		{
			Matrix mul = Matrix<T>(_matrix);
			mul = mul * (T)(-1);
			return mul;
		}

		// Operations with other Matrices

		Matrix<T> operator+ (Matrix<T> other)
		{
			size_t m = n_rows();
			size_t n = n_columns();
			if (m != other.n_rows()) throw linear_algebra_exception("Sum error: matrices have different number of rows.");
			if (n != other.n_columns()) throw linear_algebra_exception("Sum error: matrices have different number of columns.");
			Matrix<T> sum = Matrix<T>(m, n);
			for (size_t i = 0; i < m; i++)
				for (size_t j = 0; j < n; j++)
					sum[i][j] = _matrix[i][j] + other[i][j];
			return sum;
		}

		void operator+= (Matrix<T> other)
		{
			*this = *this + other;
		}

		Matrix<T> operator- (Matrix<T> other)
		{
			return (*this) + (-other);
		}

		void operator-= (Matrix<T> other)
		{
			return *this - other;
		}

		Matrix<T> operator* (Matrix<T> other)
		{
			size_t m = n_rows();
			size_t n = n_columns();
			size_t q = other.n_columns();
			if (other.isScalar())
				return (*this) * other[0][0];
			else if (isScalar())
				return other * _matrix[0][0];
			else if (n != other.n_rows()) throw linear_algebra_exception("Multiplication error: matrices do not have appropriate dimensions.");
			else
			{
				Matrix mul = Matrix<T>(m, q);
				for (size_t i = 0; i < m; i++)
					for (size_t j = 0; j < q; j++)
						for (size_t k = 0; k < n; k++)
							mul[i][j] += _matrix[i][k] * other[k][j];
				return mul;
			}
		}

		void operator*= (Matrix<T> other)
		{
			size_t m = n_rows();
			if (m != n_columns())
				throw linear_algebra_exception("Multiplication error: assigning result of multiplication between non-square matrices to self matrix.");
			if (m != other.n_rows() || m != other.n_columns())
				throw linear_algebra_exception("Multiplication error: assining result of non-square matrix multiplication to left-hand-side matrix.");
			(*this) = (*this) * other;
		}

		// Operations with alg::vectors

		vector<T> operator+ (vector<T> v)
		{
			size_t m = n_rows();
			size_t n = n_columns();
			size_t max = (m > n) ? m : n;
			size_t len = v.size();
			if (m != 1 && n != 1)
			{
				std::ostringstream error_msg("Error: sum of alg::vector with multidimensional matrix. Dimensions: ");
				error_msg << m << "x" << n;
				throw linear_algebra_exception(error_msg.str().data());
			}
			if (max != len)
			{
				std::ostringstream error_msg("Error: sum of alg::vector with unidimensional matrix of different length. Matrix length: ");
				error_msg << max << ", std::vector length: " << len;
				throw linear_algebra_exception(error_msg.str().data());
			}
			vector<T> sum = vector<T>(*this);
			for (size_t k = 0; k < max; k++)
				sum[k] += v[k];
			return sum;
		}

		void operator+= (vector<T> v)
		{
			*this = *this + v;
		}

		friend vector<T> operator+ (vector<T> v, Matrix<T> m)
		{
			return m + v;
		}

		Matrix<T> operator- (vector<T> v)
		{
			return (*this) + (-v);
		}

		void operator-= (std::vector<T> v)
		{
			*this = *this - v;
		}

		friend Matrix<T> operator- (std::vector<T> v, Matrix<T> m)
		{
			return m - v;
		}

		// If an alg::vector is multiplied by a Matrix from the right, it will be treated as a column vector.
		// The return type will still be an alg::vector with no information on whether it's supposed to be column or row.
		vector<T> operator* (vector<T> v)
		{
			Matrix vm = Matrix<T>(v); // make v into a column-vector in matrix representation.
			Matrix mul = (*this) * vm;
			return mul.toVector();
		}

		// If a Matrix is multiplied by an alg::vector from the right, the alg::vector will be treated as a row vector.
		// The return type will still be an alg::vector with no information on whether it's supposed to be column or row.
		friend vector<T> operator* (vector<T> v, Matrix<T> m)
		{
			Matrix vm = Matrix(v, false); // make v into a row-vector in matrix representation.
			Matrix mul = vm * m;
			return mul.toVector();
		}

		// Operations with scalars

		Matrix<T> operator+ (T t)
		{
			Matrix sum = Matrix<T>(*this);
			size_t m = n_rows();
			size_t n = n_columns();
			for (size_t i = 0; i < m; i++)
				for (size_t j = 0; j < n; j++)
					sum[i][j] += t;
			return sum;
		}

		friend Matrix<T> operator+ (T t, Matrix<T> M)
		{
			return M + t;
		}

		void operator+= (T t)
		{
			(*this) = (*this) + t;
		}

		Matrix<T> operator- (T t)
		{
			return (*this) + (-t);
		}

		friend Matrix<T> operator- (T t, Matrix<T> M)
		{
			return -M + t;
		}

		void operator-= (T t)
		{
			(*this) = (*this) - t;
		}

		Matrix<T> operator* (T t)
		{
			size_t m = n_rows();
			size_t n = n_columns();
			Matrix mul = Matrix<T>(_matrix);
			for (size_t i = 0; i < m; i++)
				for (size_t j = 0; j < n; j++)
					mul[i][j] *= t;
			return mul;
		}

		friend Matrix<T> operator* (T t, Matrix<T> M)
		{
			return M * t;
		}

		void operator*= (T t)
		{
			(*this) = (*this) * t;
		}

		Matrix<T> operator/ (T t)
		{
			return (*this) * (1 / t);
		}

		void operator/= (T t)
		{
			*this = (*this) / t;
		}

	};

	// Generates an identity matrix of size m. 
	template <typename T>
	Matrix<T> I(size_t m)
	{
		Matrix<T> e = Matrix<T>(m, m);
		for (size_t i = 0; i < m; i++) e[i][i] = (T)1;
		return e;
	}

}

// Overload allows for easy printing of std::vectors' contents.
template <typename T>
std::ostream& operator<< (std::ostream& output, std::vector<T> vector)
{
	size_t n = vector.size();
	std::ostringstream line;
	line << "{ ";
	for (size_t j = 0; j < n; j++)
	{
		if (j < n - 1) line << vector[j] << ", ";
		else line << vector[j] << " }";
	}
	output << line.str();
	return output;
}
