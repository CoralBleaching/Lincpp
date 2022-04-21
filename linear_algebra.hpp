#pragma once

#include "matrix.hpp"

namespace alg {

template<class T>
class Vector : private Matrix<T>
{
public:

	Vector(size_t k = 0, T val = 0)             { this->matrix_ = std::vector<T>(k, val); }
	Vector(const std::initializer_list<T> list) { this->matrix_ = std::vector<T>{ list }; }
	Vector(const std::vector<T> v)              { this->matrix_ = v; }
	Vector(const Matrix<T>& M)                  { this->matrix_ = M.getInternalStdVector(); }

	// Shape methods

	inline size_t size() const;
	inline bool isScalar() const;
	constexpr bool isVector() const;
	inline void clear();

	inline Matrix<T> row() const;
	friend Matrix<T> row(Vector<T> v);
	inline Matrix<T> col() const;
	friend Matrix<T> col(Vector<T> v);
	Matrix<T> t() const;
	inline friend Matrix<T> t(Vector<T> v);

	// Iterators and access

	inline void push_back(const T& val);
	inline void insert(Matrix<T>::Iterator it, const T& val);
	//Vector<T> concatenate(const Vector<T>&) const;
	//Matrix<T> concatenate(const Matrix<T>&) const;
	//friend Matrix<T> concatenate(const Matrix<T>&, Vector<T>&);
	Vector<T> slice(size_t start, size_t finish);

	inline typename Matrix<T>::Iterator begin() const;
	inline typename Matrix<T>::Iterator end() const;

	T& operator[](size_t i);
	T& at(size_t i);

	// Algebraic methods

	inline T norm(T) const;
	friend T norm(Vector<T> v, T val);

	// Operators

	// Unary operation

	Vector<T> operator-() const;

	// Operations with scalar

	Vector<T> operator+(T t) const;
	inline void operator+=(T t);
	friend Vector<T> operator+(T t, Vector<T> v);

	Vector<T> operator-(T t) const;
	inline void operator-=(T t);
	friend Vector<T> operator-(T t, Vector<T> v);

	Vector<T> operator*(T t) const;
	inline void operator*=(T t);
	friend Vector<T> operator*(T t, Vector<T> v);

	Vector<T> operator/(T t) const;
	inline void operator/=(T t);
	friend Vector<T> operator/(T t, Vector<T> v);

	// Operations with Vector

	T operator*(Vector<T> v) const;
	Vector<T> operator+ (const Vector<T>&) const;
	void operator+= (const Vector<T>&);
	Vector<T> operator- (const Vector<T>&) const;
	void operator-= (const Vector<T>&);

	// Operations with Matrix

	Vector<T> operator+ (const Matrix<T>&) const;
	void operator+= (const Matrix<T>& M);
	friend Vector<T> operator+(const Matrix<T>& M, Vector<T>& v);
	Vector<T> operator- (const Matrix<T>&) const;
	void operator-= (const Matrix<T>& M);
	friend Vector<T> operator-(const Matrix<T>& M, Vector<T>& v);
	Vector<T> operator*(Matrix<T> M) const;
	friend Vector<T> operator*(Matrix<T> M, Vector<T> v) { return v * M; };

};

// DEFINITIONS

template<class T>
std::ostream& operator<< (std::ostream& output, Vector<T> v) { output << toString(v); return output; }

// Constructors

// Shape methods

template<class T> inline    size_t    Vector<T>::size()		const { return Matrix<T>::size(); }
template<class T> inline    bool      Vector<T>::isScalar() const { return size() == 1; }
template<class T> constexpr bool      Vector<T>::isVector() const { return true; }
template<class T> inline    void      Vector<T>::clear()		  { Matrix<T>::clear(); }

template<class T> inline    Matrix<T> Vector<T>::row()		const { return Matrix<T>{this->matrix_, 1, size() }; }
template<class T>           Matrix<T>            row(Vector<T> v) { return v.row(); }
template<class T> inline    Matrix<T> Vector<T>::col()		const { return Matrix<T>{this->matrix_, size(), 1 }; }
template<class T>           Matrix<T>			 col(Vector<T> v) { return v.col(); }
template<class T>           Matrix<T> Vector<T>::t()		const { return col(); }
template<class T> inline    Matrix<T>	         t(Vector<T> v)	  { return v.col(); }

// Iterators and access

template<class T> inline void Vector<T>::push_back(const T& val) { this->matrix_.push_back(val); }
template<class T> void Vector<T>::insert(Matrix<T>::Iterator it, const T& val) { Matrix<T>::insert(it, val); }
template<class T> Vector<T> Vector<T>::slice(size_t start, size_t finish) 
{
	size_t n = finish - start;
	if (n > size()) throw LinearAlgebraException("Slice exceeds Vector length.");
	Vector v(n);
	std::copy_n(begin() + start, n, v.begin());
	return v;
}

template<class T> inline typename Matrix<T>::Iterator Vector<T>::begin() const { return Matrix<T>::begin(); }
template<class T> inline typename Matrix<T>::Iterator Vector<T>::end() const { return Matrix<T>::end(); }

template<class T> T& Vector<T>::operator[](size_t i) { return this->matrix_[i]; }
template<class T> T& Vector<T>::at(size_t i) { return this->matrix_.at(i); }

// Algebraic methods

template<class T> inline T Vector<T>::norm(T val = 2) const { return Matrix<T>::norm(val); }
template<class T> T norm(Vector<T> v, T val = 2) { return v.norm(val); }

// Operators

// Unary operator

template<class T>
Vector<T> Vector<T>::operator-() const
{
	Vector v(size());
	std::transform(begin(), end(), v.begin(), std::negate<T>());
	return v;
}

// Operations with scalar

template<class T>
inline Vector<T> Vector<T>::operator+(T t) const 
{
	Vector v{ *this };
	for (T& e : v) e += t; 
	return v;
}
template<class T> inline void Vector<T>::operator+=(T t) { for (T& e : this->matrix_) e += t; }
template<class T> Vector<T> operator+(T t, Vector<T> v) { return v + t; }

template<class T>
inline Vector<T> Vector<T>::operator-(T t) const
{
	Vector v{ *this };
	for (T& e : v) e -= t;
	return v;
}
template<class T> inline void Vector<T>::operator-=(T t) { for (T& e : this->matrix_) e -= t; }
template<class T> Vector<T> operator-(T t, Vector<T> v) { return v - t; }

template<class T> inline Vector<T> Vector<T>::operator*(T t) const 
{
	Vector v{ *this };
	for (T& e : v) e *= t;
	return v;
}

template<class T> inline void Vector<T>::operator*=(T t) { for (T& e : this->matrix_) e *= t; }
template<class T> Vector<T> operator*(T t, Vector<T> v) { return v * t; }

template<class T>
inline Vector<T> Vector<T>::operator/(T t) const
{
	Vector v{ *this };
	for (T& e : v) e /= t;
	return v;
}
template<class T> inline void Vector<T>::operator/=(T t) { for (T& e : this->matrix_) e /= t; }
template<class T> Vector<T> operator/(T t, Vector<T> v) { return v / t; }

// Operations with Vector

template<class T>
inline T Vector<T>::operator*(Vector<T> v) const
{
	return std::inner_product(begin(), end(), v.begin(), 0);
}

template<class T>
inline Vector<T> Vector<T>::operator+(const Vector<T>& other) const
{
	if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
	Vector v(size());
	std::transform(begin(), end(), other.begin(), v.begin(), std::plus<T>());
	return v;
}

template<class T>
inline void Vector<T>::operator+=(const Vector<T>& other)
{
	if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
	std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
}

template<class T>
inline Vector<T> Vector<T>::operator-(const Vector<T>& other) const
{
	if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
	Vector v(size());
	std::transform(begin(), end(), other.begin(), v.begin(), std::minus<T>());
	return v;
}

template<class T>
inline void Vector<T>::operator-=(const Vector<T>& other)
{
	if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
	std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
}


// Operations with Matrix

template<class T>
inline Vector<T> Vector<T>::operator+(const Matrix<T>& other) const
{
	if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
	if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
	Vector v(size());
	std::transform(begin(), end(), other.begin(), v.begin(), std::plus<T>());
	return v;
}

template<class T>
inline void Vector<T>::operator+=(const Matrix<T>& other)
{
	if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
	if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
	std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
}

template<class T> Vector<T> operator+(const Matrix<T>& M, Vector<T>& v) { return v + M; }

template<class T>
inline Vector<T> Vector<T>::operator-(const Matrix<T>& other) const
{
	if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
	if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
	Vector v(size());
	std::transform(begin(), end(), other.begin(), v.begin(), std::minus<T>());
	return v;
}

template<class T>
inline void Vector<T>::operator-=(const Matrix<T>& other)
{
	if (!other.isVector()) throw LinearAlgebraException("Can't add to Vector: Matrix isn't row or column.");
	if (size() != other.size()) throw LinearAlgebraException("Can't add Vector to Matrix of different length.");
	std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
}

template<class T> Vector<T> operator-(const Matrix<T>& M, Vector<T>& v) { return v - M; }

template<class T>
inline Vector<T> Vector<T>::operator*(Matrix<T> M) const
{
	//if (size() != M.nrows()) throw LinearAlgebraException("Can't perform dot product: Vector and Matrix have incompatible shapes.");
	Vector v;
	//for (size_t i = 0; i < M.nrows(); i++)
	//	v.push_back(std::inner_product(M.row(i).begin(), M.row(i).end(), begin(), (T)0));
	std::for_each(M.beginRow(), M.endRow(), [&](Matrix::Row& row)
		{
			v.push_back(std::inner_product(row.begin(), row.end(), begin(), 0));
		});
	return v;
}

//template<class T> Vector<T> operator*(Matrix<T> M, Vector<T> v) { return v * M; }

} // namespace::alg