#include<vector>
#include<new>
#include<sstream>
#include<ostream>
#include<cmath>

template <typename T>
std::vector<T> range(size_t);

template <typename T>
std::vector<T> range(size_t, size_t);

template <typename T>
std::ostream& operator<< (std::ostream&, std::vector<T>);

template <typename T>
std::vector<T> operator- (std::vector<T>);

template <typename T>
T operator* (std::vector<T>, std::vector<T>);

class linear_algebra_exception : public std::exception
{
	public:

	const char* message;

	linear_algebra_exception (const char * message = "")
	{
		this->message = message;
	}

	virtual const char* what()  const throw()
    {
        return this->message;
    }
};

template <class T>
class Matrix {

private:

	std::vector<std::vector<T>> matrix;

public:

	Matrix()
	{
		this->matrix = std::vector<std::vector<T>>();
	}

	Matrix(std::vector<T> user_v, bool column = true) 
	{
		if (column)
		{
			size_t m = user_v.size();
			this->matrix = std::vector<std::vector<T>>(m);
			for (int i = 0; i < m; i++) (*this)[i] = { user_v[i] };
		}
		else
			this->matrix = { user_v };
	}

	Matrix(std::vector<std::vector<T>> user_m) 
	{
		this->matrix = user_m;
	}

	Matrix(Matrix<T>* user_m)
	{
		this->matrix = std::vector<std::vector<T>>();
		size_t m = user_m->n_rows();
		for (size_t i = 0; i < m; i++) this->matrix.push_back(std::vector<T>((*user_m)[i]));
	}

	Matrix(size_t m, size_t n)
	{
		this->matrix = std::vector<std::vector<T>>(m);
		for (size_t i = 0; i < m; i++)
		{
			this->matrix[i] = std::vector<T>(n, 0);
		}
	}

	auto end()
	{
		return this->matrix.end();
	}

	T scalar()
	{
		return this->matrix[0][0];
	}

	size_t n_rows()
	{
		return this->matrix.size();
	}

	size_t n_columns()
	{
		return this->matrix[0].size();
	}

	std::vector<T> row(size_t m)
	{
		if (this->n_rows() < m) throw linear_algebra_exception("Access error: trying to access row out of index.");
		return this->matrix[m];
	}

	std::vector<T> column(size_t n)
	{
		if (this->n_columns() < n) throw linear_algebra_exception("Access error: trying to access column out of index.");
		size_t m = this->n_rows();
		std::vector<T> col(m);
		for (size_t i = 0; i < m; i++) col[i] = (*this)[i][n];
		return col;
	}

	Matrix<T> concatenate(Matrix<T> other)
	{
		Matrix M = Matrix(*this);
		size_t m = this->n_rows();
		if (m != other.n_rows()) throw linear_algebra_exception("Matrices of different number of rows.");
		for (size_t i = 0; i < m; i++)
			M[i].insert(M[i].end(), other[i].begin(), other[i].end());
		return M;
	}

	Matrix<T> slice(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end)
	{
		size_t m = row_end - row_begin;
		size_t n = column_end - column_begin;
		if (row_end > this->n_rows() || column_end > this->n_columns()) throw linear_algebra_exception("Exceeded matrix dimensions.");
		Matrix M = Matrix(m, n);
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				M[i][j] = (*this)[i + row_begin][j + column_begin];
		return M;
	}

	Matrix<T> I(size_t m)
	{
		Matrix<T> e = Matrix<T>(m, m);
		for (size_t i = 0; i < m; i++) e[i][i] = (T)1;
		return e;
	}

	Matrix<T> inverse()
	{
		size_t m = this->n_rows();
		size_t this_n = this->n_columns();
		if (m != this_n) throw linear_algebra_exception("Matrix is not square.");
		Matrix M = Matrix<T>(*this);
		M = M.concatenate(I(m));
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
		for (long long int j = m - 1; j >= 0; j--)
		{
			T mjj = 1 / M[j][j];
			for (size_t k = j; k < n; k++)
				M[j][k] *= mjj;
			for (long long int i = j - 1; i >= 0; i--)
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

	std::vector<T>& operator[] (size_t i)
	{
		return this->matrix[i];
	}

	T& operator() (size_t i, size_t j)
	{
		return this->matrix[i][j];
	}

	Matrix<T> operator+ (Matrix<T> other)
	{
		size_t m = this->n_rows();
		size_t n = this->n_columns();
		if (m != other.n_rows()) throw linear_algebra_exception("Sum error: matrices have different number of rows.");
		if (n != other.n_columns()) throw linear_algebra_exception("Sum error: matrices have different number of columns.");
		Matrix<T> sum = Matrix<T>(m, n);
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				sum[i][j] = (*this)[i][j] + other[i][j];
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

	Matrix<T> operator+ (std::vector<T> v)
	{
		size_t m = this->n_rows();
		size_t n = this->n_columns();
		size_t max = (m > n) ? m : n;
		size_t len = v.size();
		if (m != 1 && n != 1)
		{
			std::ostringstream error_msg("Error: sum of std::vector with multidimensional matrix. Dimensions: ");
			error_msg << m << "x" << n;			
			throw linear_algebra_exception(error_msg.str().data());
		}
		if (max != len)
		{
			std::ostringstream error_msg("Error: sum of std::vector with unidimensional matrix of different lenghts. Matrix length: ");
			error_msg << max << ", std::vector length: " << len;
			throw linear_algebra_exception(error_msg.str().data());
		}
		Matrix sum = Matrix<T>(*this);
		for (size_t k = 0; k < max; k++)
		{
			if (m > n)
				sum[k][0] += v[k];
			else
				sum[0][k] += v[k];
		}
		return sum;
	}

	void operator+= (std::vector<T> v)
	{
		*this = *this + v;
	}

	friend Matrix<T> operator+ (std::vector<T> v, Matrix<T> m)
	{
		return m + v;
	}

	Matrix<T> operator- (std::vector<T> v)
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

	T operator+ (T t)
	{
		if (this->isScalar()) return this->matrix[0][0] + t;
		else throw linear_algebra_exception("Error: sum of scalar with non-scalar matrix.");
	}

	T operator- (T t)
	{
		return (*this) + (-t);
	}

	void operator+= (T t)
	{
		(*this) = (*this) + t;
	}

	void operator-= (T t)
	{
		(*this) = (*this) - t;
	}

	friend T operator+ (T t, Matrix<T> M)
	{
		return M + t;
	}

	friend T operator- (T t, Matrix<T> M)
	{
		return -M + t;
	}

	friend Matrix<T> t(Matrix<T> M)
	{
		if (M.isVector())
		{
			if (M.isRow())
				return Matrix(M[0], true);
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

	friend bool isScalar(std::vector<T> v)
	{
		return v.size() == 1;
	}

	bool isScalar()
	{
		if (this->n_columns() == 1 && this->n_rows() == 1)
			return true;
		return false;
	}

	bool isColumn()
	{
		return this->n_columns() == 1;
	}

	bool isRow()
	{
		return this->n_rows() == 1;
	}

	bool isVector()
	{
		return this->isRow() || this->isColumn();
	}

	Matrix<T> operator* (Matrix<T> other)
	{
		size_t m = this->n_rows();
		size_t n = this->n_columns();
		size_t q = other.n_columns();
		if (other.isScalar())
			return (*this) * other[0][0];
		else if (this->isScalar())
			return other * (*this)[0][0];
		else if (n != other.n_rows()) throw linear_algebra_exception("Multiplication error: matrices do not have appropriate dimensions.");
		else
		{
			Matrix mul = Matrix<T>(m, q);
			for (size_t i = 0; i < m; i++)
				for (size_t j = 0; j < q; j++)
					for (size_t k = 0; k < n; k++)
						mul[i][j] += (*this)[i][k] * other[k][j];
			return mul;
		}
	}

	void operator*= (Matrix<T> other)
	{
		size_t m = this->n_rows();
		if (m != this->n_columns())
			throw linear_algebra_exception("Multiplication error: assigning result of multiplication between non-square matrices to self matrix.");
		if (m != other.n_rows() || m != other.n_columns())
			throw linear_algebra_exception("Multiplication error: assining result of non-square matrix multiplication to left-hand-side matrix.");
		(*this) = (*this) * other;
	}

	std::vector<T> operator* (std::vector<T> v)
	{
		Matrix vm = Matrix<T>(v);
		return t(*this * vm).row(0);
	}

	friend std::vector<T> operator* (std::vector<T> v, Matrix<T> m)
	{
		Matrix vm = Matrix(v, false);
		Matrix mul = vm * m ;
		return mul[0];
	}

	Matrix<T> operator* (T t)
	{
		size_t m = this->n_rows();
		size_t n = this->n_columns();
		Matrix mul = Matrix<T>(this->matrix);
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				mul[i][j] *= t;
		return mul;
	}

	Matrix<T> operator/ (T t)
	{
		return (*this) * (1 / t);
	}

	void operator*= (T t)
	{
		(*this) = (*this) * t;
	}

	void operator/= (T t)
	{
		*this = (*this) / t;
	}

	friend Matrix<T> operator* (T t, Matrix<T> matrix)
	{
		return matrix * t;
	}

	Matrix<T> operator- ()
	{
		Matrix mul = Matrix<T>(this->matrix);
		mul = mul * (T)(-1);
		return mul;
	}

	void push_back(std::vector<T> v)
	{
		if (this->matrix.size() > 0 && v.size() != this->n_columns()) 
			throw linear_algebra_exception("Error: pushing back std::vector with inappropriate dimensions.");
		this->matrix.push_back(v);
	}

	double norm()
	{
		if (this->isVector())
		{
			return (T)sqrt(((*this) * (*this))[0][0]);
		}
		throw linear_algebra_exception("Tried to take the norm of multidimensional matrix.");
	}

	static double norm(std::vector<T> v)
	{
		return sqrt(v * v);
	}

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
};

template <typename T>
Matrix<T> I(size_t m)
{
	Matrix<T> e = Matrix<T>(m, m);
	for (size_t i = 0; i < m; i++) e[i][i] = (T)1;
	return e;
}

template <typename T>
 std::vector<T> operator* (std::vector<T> v, T t)
{
	std::vector<T> mul = std::vector<T>(v);
	for (auto& x : mul)
		x *= t;
	return mul;
}

 template <typename T>
 std::vector<T> operator* (T t, std::vector<T> v)
{
	return v * t;
}

template <typename T>
T operator* (std::vector<T> lhs, std::vector<T> rhs)
{
	size_t m = lhs.size();
	if (m != rhs.size()) throw linear_algebra_exception("Error: product of two std::vector<T>'s of different lenghts.");
	T sum = 0;
	for (size_t i = 0; i < m; i++) sum += lhs[i] * rhs[i];
	return sum;
}

template <typename T>
std::vector<T> operator+ (std::vector<T> lhs, std::vector<T> rhs)
{
	size_t m = lhs.size();
	if (m != rhs.size()) throw linear_algebra_exception("Error: sum of two std::vector<T> of different lenghts.");
	std::vector<T> sum = std::vector<T>(m, 0);
	for (size_t i = 0; i < m; i++)
		sum[i] = lhs[i] + rhs[i];
	return sum;
}

template <typename T>
std::vector<T> operator- (std::vector<T> lhs, std::vector<T> rhs)
{
	return lhs + (-rhs);
}

template <typename T>
void operator+= (std::vector<T>& lhs, std::vector<T> rhs)
{
	lhs = lhs + rhs;
}

template <typename T>
void operator-= (std::vector<T>& lhs, std::vector<T> rhs)
{
	lhs = lhs - rhs;
}

template <typename T>
std::vector<T> operator+ (std::vector<T> v, T t)
{
	std::vector<T> sum = std::vector<T>(v);
	for (auto& i : sum)
		i += t;
	return sum;
}

template <typename T>
void operator+= (std::vector<T>& lhs, T& rhs)
{
	lhs = lhs + rhs;
}

template <typename T>
void operator-= (std::vector<T>& lhs, T& rhs)
{
	lhs = lhs - rhs;
}

template <typename T>
std::vector<T> operator- (std::vector<T> v, T t)
{
	return v + (-t);
}

template <typename T>
std::vector<T> operator+ (T t, std::vector<T> v)
{
	return v + t;
}

template <typename T>
std::vector<T> operator- (T t, std::vector<T> v)
{
	return v - t;
}

template <typename T>
std::vector<T> operator- (std::vector<T> v)
{
	return v * (T)(-1);
}

template <typename T>
std::vector<T> range(unsigned int n)
{
	std::vector<T> v(n);
	for (size_t i = 0; i < n; i++) v[i] = (T)i;
	return v;
}

template <typename T>
std::vector<T> range(size_t s, size_t e)
{
	std::vector<T> v(e - s);
	for (size_t i = 0; i < e - s; i++) v[i] = (T)(s + i);
	return v;
}

template <typename T>
Matrix<T> t(std::vector<T> v)
{
	return Matrix<T>(v, true);
}

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