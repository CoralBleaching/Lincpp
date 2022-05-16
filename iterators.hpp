#pragma once

namespace alg {

	template<class T> class Iterator;
	template<class T> class IteratorColumn;
	template<class T> class Row;
	template<class T> class Column;
	template<class T> class IteratorRowVector;
	template<class T> class IteratorColumnVector;

	template<class T> class Matrix;
	struct Shape { size_t m, n; } shape_;

	template<class T> 
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

	template<class T>
	class IteratorColumn
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using reference = T&;

		IteratorColumn(pointer ptr, size_t ncols) : ncols_{ ncols } { mptr = ptr; }
		IteratorColumn(Iterator<T> it, size_t ncols) : ncols_{ ncols } { mptr = it.getPtr(); }
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

	template<class T> class Vector;

	template<class T>
	class Row
	{
	public:
		Row(Iterator<T> begin = Iterator{ nullptr },
			Iterator<T> end = Iterator{ nullptr },
			Shape shape = { 0,0 }) :
			begin_{ begin }, end_{ end }, shape_{ shape } {} // automatically generate end_?

		//void operator=(Vector<T> v);

		Iterator<T> begin() const { return begin_; }
		Iterator<T> end() const { return end_; }

		inline Shape getShape() const { return shape_; }
		inline size_t ncols() const { return shape_.n; }
		inline size_t size() const { return shape_.n; }

		T& operator[](size_t i) const { return *(begin_.getPtr() + i); } //THROW
		friend std::ostream& operator<< (std::ostream& output, Row r) { output << toString(r); return output; }


		// Algebraic methods

		inline T norm(T power = 2) const;
		friend T norm(Row<T> v, T val = 2);

		// Operators

		// Unary operation

		Vector<T> operator-() const;

		// Operations with scalar

		Vector<T> operator+(T t) const;
		inline void operator+=(T t);
		friend Vector<T> operator+(T t, Row<T> v);

		Vector<T> operator-(T t) const;
		inline void operator-=(T t);
		friend Vector<T> operator-(T t, Row<T> v);

		Vector<T> operator*(T t) const;
		inline void operator*=(T t);
		friend Vector<T> operator*(T t, Row<T> v);

		Vector<T> operator/(T t) const;
		inline void operator/=(T t);
		friend Vector<T> operator/(T t, Row<T> v);

		// Operations with Vector

		T operator*(Vector<T> v) const;
		friend T operator*(Vector<T> v, Row<T> r) { return r * v; }
		Vector<T> operator+ (const Vector<T>&) const;
		void operator+= (const Vector<T>&);
		friend Vector<T> operator+(Vector<T> v, Row<T> r) { return r + v; }
		Vector<T> operator- (const Vector<T>&) const;
		void operator-= (const Vector<T>&);
		friend Vector<T> operator-(Vector<T> v, Row<T> r) { return r - v; }

		// Operations with Row

		T operator*(Row<T> v) const;
		Vector<T> operator+ (const Row<T>&) const;
		void operator+= (const Row<T>&);
		Vector<T> operator- (const Row<T>&) const;
		void operator-= (const Row<T>&);

	protected:
		Iterator<T> begin_, end_;
		Shape shape_;
	};

	template<class T>
	class Column
	{
	public:
		Column(Iterator<T> begin = Iterator{ nullptr },
			Iterator<T> end = Iterator{ nullptr },
			Shape shape = { 0,0 }) :
			begin_{ IteratorColumn<T>{begin, shape.n} }, end_{ IteratorColumn<T>{end, shape.n} }, shape_{ shape } {}

		IteratorColumn<T> begin() const { return begin_; }
		IteratorColumn<T> end() const { return end_; }

		inline Shape getShape() const { return shape_; }
		inline size_t nrows() const { return shape_.m; }
		inline size_t size() const { return shape_.m; }

		T& operator[](size_t i) { return *(begin_ + i); } //THROW
		friend std::ostream& operator<< (std::ostream& output, Column c) { output << toString(c); return output; }

		// Algebraic methods

		inline T norm(T power = 2) const;
		friend T norm(Column<T> v, T val = 2);

		// Operators

		// Unary operation

		Vector<T> operator-() const;

		// Operations with scalar

		Vector<T> operator+(T t) const;
		inline void operator+=(T t);
		friend Vector<T> operator+(T t, Column<T> v);

		Vector<T> operator-(T t) const;
		inline void operator-=(T t);
		friend Vector<T> operator-(T t, Column<T> v);

		Vector<T> operator*(T t) const;
		inline void operator*=(T t);
		friend Vector<T> operator*(T t, Column<T> v);

		Vector<T> operator/(T t) const;
		inline void operator/=(T t);
		friend Vector<T> operator/(T t, Column<T> v);

		// Operations with Vector

		T operator*(Vector<T> v) const;
		friend T operator*(Vector<T> v, Column<T> c) { return c * v; }
		Vector<T> operator+ (const Vector<T>&) const;
		void operator+= (const Vector<T>&);
		friend Vector<T> operator+(Vector<T> v, Column<T> c) { return c + v; }
		Vector<T> operator- (const Vector<T>&) const;
		void operator-= (const Vector<T>&);
		friend Vector<T> operator-(Vector<T> v, Column<T> c) { return c - v; }

		// Operations with Column

		T operator*(Column<T> v) const;
		Vector<T> operator+ (const Column<T>&) const;
		void operator+= (const Column<T>&);
		Vector<T> operator- (const Column<T>&) const;
		void operator-= (const Column<T>&);

	private:
		IteratorColumn<T> begin_, end_;
		Shape shape_;
	};

	template<class T>
	class IteratorRowVector
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = Row<T>;
		using pointer = Iterator<T>;
		using reference = Row<T>&;

		IteratorRowVector(Iterator<T> it = Iterator<T>{ nullptr }, Shape shape = { 0,0 }) : it_{ it }, shape_{ shape } {}
		IteratorRowVector(const IteratorRowVector& rawIterator) = default;
		IteratorRowVector& operator=(const IteratorRowVector& rawIterator) = default;

		Row<T>& operator*() { Row<T> r{ it_, it_ + ncols(), shape_ }; return std::move(r); }

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

		Iterator<T> getIt() { return it_; }
		size_t nrows() { return shape_.m; }
		size_t ncols() { return shape_.n; }

	private:
		Iterator<T> it_;
		Shape shape_;
	};

	template<class T>
	class IteratorColumnVector
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;

		IteratorColumnVector(const Matrix<T>* ptr = nullptr, size_t index = 0, Shape shape = { 0,0 }) : mptr{ ptr }, index_{ index }, shape_{ shape } {}

		Column<T>& operator*()
		{
			Column<T> c = mptr->col(index_);
			return std::move(c);
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

	// Row

	// Algebraic methods

	template<class T> inline T Row<T>::norm(T power = 2) const
	{
		return sqrt(std::accumulate<Iterator<T>, T>(begin(), end(), 0, [=](T& accum, T& next)
			{
				return accum + pow(next, power);
			}));
	}
	template<class T> T norm(Row<T> v, T val = 2) { return v.norm(val); }

	// Operators

	// Unary operator

	template<class T>
	Vector<T> Row<T>::operator-() const
	{
		Vector<T> v(size());
		std::transform(begin(), end(), v.begin(), std::negate<T>());
		return v;
	}

	// Operations with scalar

	template<class T>
	inline Vector<T> Row<T>::operator+(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e += t;
		return v;
	}
	template<class T> inline void Row<T>::operator+=(T t) 
	{
		std::for_each(begin(), end(), [&](T& e) { e += t; });
	}
	template<class T> Vector<T> operator+(T t, Row<T> v) { return v + t; }

	template<class T>
	inline Vector<T> Row<T>::operator-(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e -= t;
		return v;
	}
	template<class T> inline void Row<T>::operator-=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e -= t; });
	}
	template<class T> Vector<T> operator-(T t, Row<T> v) { return v - t; }

	template<class T> inline Vector<T> Row<T>::operator*(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e *= t;
		return v;
	}

	template<class T> inline void Row<T>::operator*=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e *= t; });
	}
	template<class T> Vector<T> operator*(T t, Row<T> v) { return v * t; }

	template<class T>
	inline Vector<T> Row<T>::operator/(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e /= t;
		return v;
	}
	template<class T> inline void Row<T>::operator/=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e /= t; });
	}
	template<class T> Vector<T> operator/(T t, Row<T> v) { return v / t; }

	// Operations with Vector

	template<class T>
	inline T Row<T>::operator*(Vector<T> v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0);
	}

	template<class T>
	inline Vector<T> Row<T>::operator+(const Vector<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<T>());
		return v;
	}

	template<class T>
	inline void Row<T>::operator+=(const Vector<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
	}

	template<class T>
	inline Vector<T> Row<T>::operator-(const Vector<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<T>());
		return v;
	}

	template<class T>
	inline void Row<T>::operator-=(const Vector<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
	}


	// Operations with Row

	template<class T>
	inline T Row<T>::operator*(Row<T> v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0);
	}

	template<class T>
	inline Vector<T> Row<T>::operator+(const Row<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<T>());
		return v;
	}

	template<class T>
	inline void Row<T>::operator+=(const Row<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
	}

	template<class T>
	inline Vector<T> Row<T>::operator-(const Row<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<T>());
		return v;
	}

	template<class T>
	inline void Row<T>::operator-=(const Row<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
	}

	// Column

	// Algebraic methods

	template<class T> inline T Column<T>::norm(T power = 2) const
	{
		return sqrt(std::accumulate<Iterator<T>, T>(begin(), end(), 0, [=](T& accum, T& next)
			{
				return accum + pow(next, power);
			}));
	}
	template<class T> T norm(Column<T> v, T val = 2) { return v.norm(val); }

	// Operators

	// Unary operator

	template<class T>
	Vector<T> Column<T>::operator-() const
	{
		Vector<T> v(size());
		std::transform(begin(), end(), v.begin(), std::negate<T>());
		return v;
	}

	// Operations with scalar

	template<class T>
	inline Vector<T> Column<T>::operator+(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e += t;
		return v;
	}
	template<class T> inline void Column<T>::operator+=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e += t; });
	}
	template<class T> Vector<T> operator+(T t, Column<T> v) { return v + t; }

	template<class T>
	inline Vector<T> Column<T>::operator-(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e -= t;
		return v;
	}
	template<class T> inline void Column<T>::operator-=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e -= t; });
	}
	template<class T> Vector<T> operator-(T t, Column<T> v) { return v - t; }

	template<class T> inline Vector<T> Column<T>::operator*(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e *= t;
		return v;
	}

	template<class T> inline void Column<T>::operator*=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e *= t; });
	}
	template<class T> Vector<T> operator*(T t, Column<T> v) { return v * t; }

	template<class T>
	inline Vector<T> Column<T>::operator/(T t) const
	{
		Vector<T> v{ begin(), end() };
		for (T& e : v) e /= t;
		return v;
	}
	template<class T> inline void Column<T>::operator/=(T t)
	{
		std::for_each(begin(), end(), [&](T& e) { e /= t; });
	}
	template<class T> Vector<T> operator/(T t, Column<T> v) { return v / t; }

	// Operations with Vector

	template<class T>
	inline T Column<T>::operator*(Vector<T> v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0);
	}

	template<class T>
	inline Vector<T> Column<T>::operator+(const Vector<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<T>());
		return v;
	}

	template<class T>
	inline void Column<T>::operator+=(const Vector<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
	}

	template<class T>
	inline Vector<T> Column<T>::operator-(const Vector<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<T>());
		return v;
	}

	template<class T>
	inline void Column<T>::operator-=(const Vector<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
	}


	// Operations with Row

	template<class T>
	inline T Column<T>::operator*(Column<T> v) const
	{
		return std::inner_product(begin(), end(), v.begin(), 0);
	}

	template<class T>
	inline Vector<T> Column<T>::operator+(const Column<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::plus<T>());
		return v;
	}

	template<class T>
	inline void Column<T>::operator+=(const Column<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't add Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
	}

	template<class T>
	inline Vector<T> Column<T>::operator-(const Column<T>& other) const
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		Vector<T> v(size());
		std::transform(begin(), end(), other.begin(), v.begin(), std::minus<T>());
		return v;
	}

	template<class T>
	inline void Column<T>::operator-=(const Column<T>& other)
	{
		if (size() != other.size()) throw LinearAlgebraException("Can't subtract Vectors of different lengths.");
		std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
	}

} // namespace::alg
