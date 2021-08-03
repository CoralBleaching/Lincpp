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

-  Direct multiplication and type handling (and other binary operations)

````C++
vector<double> v = { 5, 6 };
vector<double> mult = A * v;
std::cout << mult;
````
Output:
````
{ 17, 39 }
````

- Unary operations such as inverse, determinant, norm etc.
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
For a full example of usage in an application, check [this project](https://github.com/CoralBleaching/Non-linear-optimization) on non-linear optimization. It makes use of this header on its core routines.

## Description
The module was design to provide natural syntax for mathematical operations on vectors and matrices, much like the Numpy package for the Python language or MATLAB. For instance, this means that the dot product between two vectors $v$ and $u$ should be written by simple use of the `*` operator for ordinary multiplication:
````C++
double w = v * u;
````
Of course, as you can see from the above example, the results of operations should have reasonable types.
> ##### Note
> Therefore, the adopted _philosophy of operations_ is that operations that can return a "lesser" type _generally_ do return a "lesser" type. 
>
> For instance, as you can see from the examples above, multiplication between matrices and vectors return vector types (`alg::vector<T>`), and the dot product returns a scalar type (of type `T`). 

**The header defines its own namespace, `alg`.**

Vector properties and behavior was encapsulated in a template class called `vector`. Matrix properties and behavior was encapsulated in a template class called `Matrix`.
> ##### Note:
> All methods return new instances, leaving the original objects intact.

Vector do not distinguish between column or row types. When couple with other objects like matrices and other vectors, the appropriate behavior is inferred. For example, both `A*v` and `v*A`, where `A` is a matrix and `v` is a vector, give the same result, as `v` is automatically transposed accordingly.

Also,
> ##### Note:
> There is no safeguard against rows having different lengths within a single matrix.

### The `vector` class


### The `Matrix` class

`Matrix<T>` defines an object representing a matrix which behaves like a `vector<vector<T>>` object, whose elements can be accessed via brackets (`matrix[i][j]` or `matrix[i]`) or parentheses (`matrix(i, j)`).

#### Constructors
| | |
|--|--|
|`Matrix()`| Default constructor |
|`Matrix(std::vector<T>, bool)`| Creates  a matrix object encapsulating a vector. If the second argument is true, the shape of the matrix is a column matrix, otherwise row (true by default).|
|`Matrix(Matrix<T>&)`| Copy constructor |
|`Matrix(unsigned int, unsigned int)`| Allocates an _m_ x _n_ matrix, with entries initialized to zero.|
|`Matrix(std::initializer_list<std::initializer_list<T>>)` | Allows initialization via initializer lists.|
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
|`isColumn` | Returns true if the matrix has a single column, false otherwise |
|`isScalar()` | Returns true if the matrix is 1x1, false otherwise |
|`isRow()` | Returns true if the matrix has a single row, false otherwise |
|`isVector()` | Returns true if the matrix has either a single column or a single row |
|`norm()` | If the matrix is either row or column, returns its norm |
|`norm(vector<T> v)` |Sobrecarga do método anterior|
|`n_rows()`| The number of rows in the matrix |
| `n_columns()` | The number of columns in the matrix |
| `push_back(vector<T>)` | Appends a new line to the bottom of the matrix |
| `row(unsigned int)` | Returns a vector containing the specified row |
|`slice(unsigned int, unsigned int, unsigned int, unsigned int)` | Returns a submatrix of this matrix. Arguments are: upper row index, bottom row index, leftmost column index and rightmost column index. Bottom and rightmost indexes are exclusive|
|`t(Matrix<T>)` | Returns the transpose of a Matrix |
|`t(vector<T>)` | Creates a column matrix containing the elements of a vector. It's the same as `Matrix<T>(v, true)` |
|`toVector()`| If the matrix has either a single column or a single row, returns a vector containing its elements (erasing information about orientation) |

#### Operators
Matematical
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

Access
+ `[]` Returns a vector object, the row at the respective index. (From _0_ to _m - 1_)
+ `()` Same as `[]`. If given two arguments, behaves the same as `[][]`. Example: `matrix(i, j) == matrix[i][j]`.

Miscellaneous
+ `<<` Appropriately sends matrix objects to the designated stream.
