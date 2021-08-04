# Linear algebra in C++

## Introduction
The `linear_algebra.hh` header implements functionality and a numerical interface for linear algebra objects and constructs. Here is a short list of some of the things you can do with the functionality defined in this header:

- Easily declare and _print_ matrix and vector objects.

````C++
#include "linear_algebra.hh"
#include<iostream>
using namespace::alg;
int main() {
	Matrix<double> A = { {1, 2}, {3, 4} };
	std::cout << A;
}
````

Output:
	
````
{{ 1, 2 },
 { 3, 4 }}
````

-  Direct multiplication (and other binary operations) and automatic type handling

````C++
vector<double> v = { 5, 6 };
vector<double> mult = A * v;
std::cout << mult;
````
Output:
````
{ 17, 39 }
````

- Unary operations and routines such as inverse, determinant, norm etc.
````C++
std::cout << A[0].norm() << "\n";
std::cout << A.inverse() << "\n";
````
Output:
````
2.23607
{{ -2, 1 },
 { 1.5, -0.5 }}
````

- Scalar operations.
````C++
std::cout << 2*A + 1 << "\n";
````
Output:
````
{{ 3, 5 },
 { 7, 9 }}
```` 
#### Example of use
For a full example of application, check [this project](https://github.com/CoralBleaching/Non-linear-optimization) on non-linear optimization. It makes use of this header on its core routines.

Below is one function defined in that project. It takes a functional object `G` that represents a vector function, for instance, the gradient of a multivariable function, _**G**(**x**)_ = ∇ _f(**x**)_, and returns another function that computes its Jacobian at any specific point via a central difference formula:
> _H<sub>ij</sub>_ = [***G***(***x*** + _e·I<sub>j</sub>_) -  ***G***(***x*** - _e·I<sub>j</sub>_)]<sub>_i_</sub> / (_4e_) +  [***G***(***x*** + _e·I<sub>i</sub>_) -  ***G***(***x*** - _e·I<sub>i</sub>_)]<sub>_j_</sub> / (_4e_),

where _H_ = [ _H<sub>ij</sub>_ ] is the resulting Jacobian matrix, _I_ is the identity matrix, and _e_ is the desired precision. (Here we call the Jacobian "H" instead of "J" because in this project the function is tasked to find the Hessian of a function _f(**x**)_ from its gradient, or its second derivative from its first derivative, in case of a _x_ being a scalar.)

````C++
matrix_function jacob(std:function<vector<double>(vector<double>)> G, int n, double e = 1e-4)
{
    return [G, n, e](vecd x) {
        Matrix<double> H(n, n);
        Matrix<double> M = e * I<double>(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                H[i][j] = (G(x + M(j)) - G(x - M(j)))(i)  / (4 * e) + (G(x + M(i)) - G(x - M(i)))(j) / (4 * e);
        H = 0.5 * (H + t(H));
        return H;
    };
} // the return type 'matrix_function' is std::function<Matrix<double>(vector<double>)>
````

## Description
The module was designed to provide a natural syntax for mathematical operations on vectors and matrices, much like the Numpy package for the Python language or MATLAB. For instance, this means that the dot product between two vectors ***v*** and ***u*** should be written with simple use of the `*` operator for ordinary multiplication:
````C++
vector<double> v = { v1, ..., vn };
vector<double> u = { u1, ..., un };
double w = v * u;
````
Of course, as you can see from the above example, the results of operations should have reasonable types.
> ##### Note
> Therefore, the adopted _philosophy of operations_ is that operations that can return a "lesser" type _generally_ do return a "lesser" type. 
>
> For instance, as you can see from the examples above, multiplication between matrices and vectors return vector types (`alg::vector<T>`), and the dot product returns a scalar type (of type `T`). 
>
> Nevertheless, there are instances when you don't want the result of a dot product to be reduced to an object of type `vector<T>`, or the result might still be a matrix, for example: ( _v<sub>1</sub>_, ..., _v<sub>m</sub>_ )<sup>_T_</sup> · ( _u<sub>1</sub>_, ..., _u<sub>m</sub>_ ). For this reason, operations between two objects of type `Matrix<T>` still always return a new object of type `Matrix<T>`.

**The header defines its own namespace, `alg`.**

Vector properties and behavior were encapsulated in a template class called `vector`. Matrix properties and behavior were encapsulated in a template class called `Matrix`.
> ##### Note:
> All methods return new instances, leaving the original objects intact.

Vectors do not distinguish between column or row types. When coupled with other objects like matrices and other vectors, the appropriate behavior is inferred. For example, both `A*v` and `v*A`, where `A` is a matrix and `v` is a vector, give the same result, as `v` is automatically transposed accordingly.

Also,
> ##### Note:
> There is no safeguard against rows having different lengths within a single matrix.

### The `vector` class

`alg::vector<T>` defines an object representing a mathematical vector which behaves in C++ like an `std::vector` with augmented operators. Its elements can be access via brackets or parentheses (`vector[i]` or `vector(i)`).

#### Constructors
|||
|-|-|
| `vector()` | Default constructor |
| `vector(vector<T>)` | Copy constructor |
| `vector(std::vector<T>)` | Casting constructor |
| `vector(Matrix<T>)` | Works if matrix is either row or column, but that information is lost |
| `vector(unsigned int, T)` | Allocates a vector of specified length and initializes all entries with the same value. Default value: 0 |
| `vector(std::initializer_list<T>` | Allows initialization via initializer lists.|

#### Methods
|||
|-|-|
|`back()`| Returns the last element by reference |
|`begin()` | Returns an iterator pointing to the first element |
|`concatenate(vector<T>)`| Concatenates the elements of another vector to the right, returning a new, longer vector |
|`end()` | Returns an iterator pointing to just past the last element |
|`insert(Iterator, Iterator, Iterator)`| Analogous to the `insert` method of `std::vector` |
|`isScalar()`| Returns true if the number of elements in vector is 1, false otherwise |
|`norm()` | The Euclidean norm of the vector |
|`pop_back()`| Removes the last element |
|`push_back(T)`| Appends a new element after the last element |
|`size()`| Returns the length of the vector |
|`slice(unsigned int, unsigned int)` | Returns a vector comprised of a sublist of elements within a range of indexes of the original vector |

#### Functions on `vector<T>`
|||
|-|-
|`norm(vector<T>)` | Returns the Euclidean norm of the given vector |
|`t(vector<T>)` | Creates a column matrix containing the elements of a vector. It's the same as `Matrix<T>(v, true)` |

#### Operators
Mathematical 
+  `+` commutative. When used with:
	+ scalars
		+  Shifts every element by the scalar.  Return type: `vector<T>`
	+ vectors
		+ Sums only if operands are of same length. Return type: `vector<T>`
+ `-` if
	+ unary:
		+ scales each element by -1.
	+ binary:
		+ same behavior as `+`.
+ `*` commutative. When used with:
	+ scalars
		+ Scales every element by the scalar. Return type: `vector<T>`
	+ vectors	
		+ If operands of same length, does the dot product and returns the result. Return type: `T`
		> It is commutative because it always treats the left-hand side operand as a row-vector and the right-hand side operand as a column vector.

+ `/` non-commutative. Can only be used with:
	+ scalars
		+ same behavior as `*`.
+ `+=` `-=` `*=` `/=` also avaliable. 
> **For operations with matrices, please refer to the section on operators of the `Matrix` class.**

Access
+ `[]` Returns a `T` object (by reference), the entry at the respective index. (From _0_ to _n - 1_).
+ `()` Same as `[]`.

Assignment  and comparison
+ `=` Appropriately copies the contents of the right-hand side vector to the left-hand side one.
+ `==` Checks for element by element equality.
+ `!=` Checks whether some elements with the same index differ between two vectors of same length.

Miscellaneous
+ `<<` Appropriately sends vector objects to the designated stream.

### The `Matrix` class

`alg::Matrix<T>` defines an object representing a matrix which behaves like a `vector<vector<T>>` object, whose elements can be accessed via brackets (`matrix[i][j]` or `matrix[i]`) or parentheses (`matrix(i, j)`).

#### Constructors
| | |
|--|--|
|`Matrix()`| Default constructor |
|`Matrix(std::vector<T>, bool)`| Creates  a matrix object encapsulating a vector. If the second argument is true, the shape of the matrix is a column matrix, otherwise row (true by default).|
|`Matrix(Matrix<T>&)`| Copy constructor |
|`Matrix(unsigned int, unsigned int)`| Allocates an _m_ x _n_ matrix, with entries initialized to zero.|
|`Matrix(std::initializer_list<std::initializer_list<T>>)` | Allows initialization via initializer lists|
|||

#### Methods
| | |
|--|--|
|`begin()` | Returns an iterator pointing to the first element |
|`column(unsigned int)`| Returns a Matrix object comprised of the specified column |
|`determinant()` | Returns the determinant of this matrix if it's a square matrix |
|`end()` | Returns an iterator pointing to just past the last element |
|`concatenate(Matrix<T>)` | Concatenates two matrices horizontally |
|`inverse()` | Returns the inverse of this matrix of this matrix.  |
|`isColumn()` | Returns true if the matrix has a single column, false otherwise |
|`isScalar()` | Returns true if the matrix is 1x1, false otherwise |
|`isRow()` | Returns true if the matrix has a single row, false otherwise |
|`isVector()` | Returns true if the matrix has either a single column or a single row |
|`norm()` | If the matrix is either row or column, returns its norm |
|`n_rows()`| The number of rows in the matrix |
| `n_columns()` | The number of columns in the matrix |
| `push_back(vector<T>)` | Appends a new line to the bottom of the matrix |
| `row(unsigned int)` | Returns a vector containing the specified row |
|`slice(unsigned int, unsigned int, unsigned int, unsigned int)` | Returns a submatrix of this matrix. Arguments are: upper row index, bottom row index, leftmost column index and rightmost column index. Bottom and rightmost indexes are exclusive|
|`toVector()`| If the matrix has either a single column or a single row, returns a vector containing its elements (erasing information about orientation) |

#### Functions on `Matrix<T>`
|||
|-|-|
|`norm(Matrix<T> m)` | Gives the norm of `m` if it's a column or row-matrix |
|`t(Matrix<T>)` | Returns the transpose of a Matrix |

#### Operators
Mathematical
+  `+` commutative. When used with:
	+ scalars
		+  Shifts every element by the scalar.  Return type: `Matrix<T>`
	+ vectors
		+ Sums only if matrix is row/column of same length as the vector. Return type: `vector<T>`
	+ matrices
		+ Sums both matrices if they are of same shape.  Return type: `Matrix<T>`
+ `-` if
	+ unary:
		+ scales each element by -1.
	+ binary:
		+ same behavior as `+`.
+ `*` commutative. When used with:
	+ scalars
		+ Scales every element by the scalar. Return type: `Matrix<T>`
	+ vectors	
		+ If operands are of appropriate shapes, does matrix multiplication and returns the resulting vector (of same length as the one given in the argument). Return type: `vector<T>`
	+ matrices
		+ Does matrix multiplication and returns a new matrix of appropriate shape.  Return type: `Matrix<T>`
+ `/` non-commutative. Can only be used with:
	+ scalars
		+ same behavior as `*`.
+ `+=` `-=` `*=` `/=` also avaliable. 

Access
+ `[]` Returns a vector object, the row at the respective index. (From _0_ to _m - 1_)
+ `()` Same as `[]`. If given two arguments, behaves the same as `[][]`. Example: `matrix(i, j) == matrix[i][j]`.

Assignment
+ `=` Appropriately copies the contents of the right-hand side matrix to the left-hand side one.

Miscellaneous
+ `<<` Appropriately sends matrix objects to the designated stream.
