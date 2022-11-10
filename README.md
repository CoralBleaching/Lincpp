# Linear algebra in C++

## Summary
- [Introduction](#introduction)
	- [Example of use](#example-of-use)
- [Description](#description)
- [Vector class](#the-vector-class)
	- [Constructors](#constructors)
	- [Methods](#methods)
	- [Functions on Vector](#functions-on-algvector)
	- [Operators](#operators)
- [Matrix class](#functions-on-matrixt)
	- [Constructors](#constructors-1)
	- [Methods](#methods-1)
	- [Functions on Vector](#functions-on-algmatrix)
	- [Operators](#operators-1)
- [Convenience functions](#convenience-functions)
- [Internal types (Row and Column)](#internal-types-algrow-and-algcolumn)

## Introduction
The `linear_algebra.hpp` header implements functionality and a numerical interface for linear algebra objects and constructs. Here is a short list of some of the things you can do with the functionality defined in this header:

- Easily declare and _print_ matrix and vector objects.

````C++
#include "linear_algebra.hpp"
#include<iostream>
using alg::Matrix;

int main() {
	Matrix A = { {1, 2}, {3, 4} };
	std::cout << A;
}
````

Output:
	
````

{{ 1 2 }
 { 3 4 }}

````

-  Direct multiplication (and other binary operations) and automatic type handling

````C++
using alg::Vector;

Vector v{ { 5, 6 } };
Vector mult = A * v;
std::cout << mult;
````
Output:
````
{ 17 39 }
````

- Unary operations and routines such as inverse, determinant, norm etc.
````C++
std::cout << A[0].norm() << "\n";
std::cout << A.inverse();
````
Output:
````
2.23607

{{ -2 1 }
 { 1.5 -0.5 }}

````

- Scalar operations.
````C++
std::cout << 2*A + 1;
````
Output:
````

{{ 3 5 }
 { 7 9 }}

```` 
#### Example of use
For a full example of application, check [this project](https://github.com/CoralBleaching/Non-linear-optimization) on non-linear optimization. It makes use of this header on its core routines.

Below is one function defined in that project. It takes a functional object `G` that represents a vector function, for instance, the gradient of a multivariable function, $G(x) = \nabla f(x)$, and returns another function that computes its Jacobian at any specific point via a central difference formula:

$$
H_{ij} = \frac{[G(\vec x + eI_j) - G(\vec x - eI_j)]_i + [G(\vec x + eI_i) - G(\vec x - eI_i)]_j}{4e},
$$

where $H = (H_{ij})$ is the resulting Jacobian matrix, $I$ is the identity matrix, and $e$ is the desired precision. Here we call the Jacobian $H$ instead of $J$ because in this project the function is tasked to find the Hessian of a function $f(x)$ from its gradient (or its second derivative from its first derivative, in case of a $x$ being a scalar).

````C++
#include <functional>

std::function<Matrix(Vector)> jacobian(std::function<Vector(Vector)> G, int n, double e = 1e-4)
{
    return [G, n, e](Vector x) {
        Matrix H(n, n);
        Matrix M = e * I(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                H[i][j] = (G(x + M.row(j)) - G(x - M.row(j))).at(i) / (4 * e) + (G(x + M.row(i)) - G(x - M.row(i))).at(j) / (4 * e);
        H = 0.5 * (H + t(H));
        return H;
    };
}
````

## Description
The module was designed to provide a natural syntax for mathematical operations on vectors and matrices, much like the Numpy package for the Python language or MATLAB. For instance, this means that the dot product between two vectors $\vec v$ and $\vec u$ should be written with simple use of the `*` operator for ordinary multiplication:
````C++
Vector v{ v1, ..., vn };    // Some different ways to initialize a
Vector u = { u1, ..., un }; // Vector with specific values.
double w = v * u;
````
Of course, as you can see from the above example, the results of operations have reasonable types.
> ##### Note
> Therefore, the adopted _philosophy of operations_ is that operations that can return a "lesser" type _generally_ do return a "lesser" type. 
>
> For instance, as you can see from the examples above, multiplication between matrices and vectors return Vector types (`alg::Vector`), and the dot product returns a scalar type (of type `double`). 
>
> Nevertheless, there are instances when you don't want the result of a dot product to be reduced to an object of type `alg::Vector`, or the result might still be a matrix, for example: $(v_1, \dots, v_m)^T\cdot(u_1,\dots,u_m)$. For this reason, operations between two objects of type `alg::Matrix` still always return a new object of type `alg::Matrix`. Operations of the type `A.col(j) * A.row(i)` (where `A` is a square `alg::Matrix`) also yield an `alg::Matrix` return type.

**The header defines its own namespace, `alg`.**

Vector properties and behavior were encapsulated in a class called `Vector`. Matrix properties and behavior were encapsulated in a class called `Matrix`.
> ##### Note:
> Methods and operators return new copies of the objects. Few methods act in-place, mostly insertion and access methods. For instance, the `reshape` method returns a copy, but the `setShape` method acts in-place.

Vectors do not distinguish between column or row types. When coupled with other objects like matrices and other vectors, the appropriate behavior is inferred. For example, both `A*v` and `v*A`, where `A` is a matrix and `v` is a vector, give the expected results, as `v` is automatically transposed accordingly.

## The `Vector` class

`alg::Vector` defines an object representing a mathematical vector which behaves in C++ like an `std::vector` with augmented operators. Its elements can be accessed via brackets or `at` method (`vector[i]` or `vector.at(i)`).

### Constructors
|||
|-|-|
| `alg::Vector(size_t, double)` | Default constructor (defaults: `0`, `0.`). Explicit. Allocates a vector of `size_t` length and initializes all entries with the same `double` value. |
| `alg::Vector(std::vector<double>)` | Conversion constructor for `std::vector`. Explicit. |
| `alg::Vector(alg::Matrix)` | Flattens the Matrix, losing shape information. Explicit. |
| `alg::Vector(alg::Iterator, alg::Iterator)` | Initializes a new `alg::Vector` with values from the range defined by the iterators. There's a variant for the built-in `alg::IteratorColumn`. |
| `alg::Vector(std::initializer_list<double>)` | Allows initialization via initializer lists.|

### Methods
|||
|-|-|
|`at(size_t)` | Returns a `const` reference to the element at given position. Will throw if position is invalid. |
|`begin()` | Returns an iterator pointing to the first element |
|`clear()` | Clears the contents of the object. Size is set to `0`. |
|`col()` | Returns a column `alg::Matrix` with the same values. |
|`concatenate(alg::Vector)`| Concatenates the elements of another vector to the right, returning a new, longer vector |
|`end()` | Returns an iterator pointing to just past the last element |
|`getInternalStdVector()`| Returns a reference to the underlying data container. |
|`insert(Iterator, double)`| Inserts a value at the position specified by the iterator. |
|`insert(Iterator, alg::Vector)`| Inserts the entire contents of a Vector at the position specified by the iterator. |
|`isScalar()`| Returns true if the number of elements in vector is 1, false otherwise. |
|`norm(double)` | The $p$-norm of the vector. Argument is the $p$ parameter. Default is `2.`. |
|`push_back(double)`| Appends a new element after the last element. Note that this changes many mathematical properties of the object. |
|`row()`| Returns a row `alg::Matrix` with the same values. Same as `t()`. |
|`size()`| Returns the length of the vector. |
|`slice(size_t, size_t)` | Returns a vector comprised of a sublist of elements within a range of indexes of the original vector. Will throw if limit exceeds the vector's lenght. |
| `t()` | Returns a row `alg::Matrix` with the same values. Same as `row()`. |
| `to_string()` | Returns a representation of the Vector as `std::string`. |

### Functions on `alg::Vector`
|||
|-|-|
|`alg::col(alg::Vector)`| Returns a column `alg::Matrix` with the same values. |
|`alg::norm(alg::Vector, double)` | The $p$-norm of the vector. Second argument is the $p$ parameter. Default is `2.`. |
|`alg::row(alg::Vector)`| Returns a row `alg::Matrix` with the same values. Same as `t(vector)`.
|`alg::t(alg::Vector)` | Returns a row `alg::Matrix` with the same values. Same as `t(vector)` or `vector.t()`. |

### Operators
Mathematical 
+  `+` commutative. When used with:
	+ scalars
		+  Shifts every element by the scalar.  Return type: `alg::Vector`
	+ vectors
		+ Sums only if operands are of same length. Return type: `alg::Vector`
	+ matrices
		+ Sums only if Matrix is row or column of same length. Return type: `alg::Vector`
+ `-` if
	+ unary:
		+ scales each element by -1.
	+ binary:
		+ same behavior as `+`.
+ `*` commutative with vectors. Non-commutative with matrices. When used with:
	+ scalars
		+ Scales every element by the scalar. Return type: `alg::Vector`
	+ vectors	
		+ If operands of same length, does the dot product and returns the result. Return type: `double`
		> It is commutative because it always treats the left-hand side operand as a row-vector and the right-hand side operand as a column vector.
	+ matrices
		+ Treats the vector as column and multiplies. May throw if shapes misalign. Return type: `alg::Vector`
+ `/` non-commutative. Can only be used with:
	+ scalars
		+ same behavior as `*`.
+ `+=` `-=` `*=` `/=` also avaliable. 

Access
+ `[]` Returns a `double` by reference, the entry at the respective index.

Other operators
+ `=` If lenghts are the same, appropriately copies the contents of the right-hand side to the left-hand side. If used with a matrix, the matrix must be either row or column of same length.
+ `<<` Appropriately sends vector objects to the designated stream.

## The `Matrix` class

`alg::Matrix` defines an object representing a matrix which behaves like an array of vector objects although it's implement with contiguous data. It uses C-style array syntax, or "linewise". Its elements can be accessed via parentheses (`matrix(i, j)`) or brackets (`matrix[i][j]`). The parentheses method is the faster one. `matrix[i]` will return a `Row` object that can be further indexed to extract an element at position `j`.

### Shape

`alg::Shape` is a struct keeping a pair of fields `m`, `n` of the number of rows and columns, respectively, of a matrix. It's shorthand used for passing this information packaged together. Like any simple object it can be instantiated with brackets (initializer lists), e.g. `{2, 3}`.

It also features a comparison operator `==` for testing equality of shapes.

### Constructors
| | |
|--|--|
|`alg::Matrix(size_t, alg::Shape)`| Default constructor. Builds a matrix with `size_t` elements and given `Shape`. Defaults are `0` and `{ 0, 0 }`, respectively. |
|`alg::Matrix(size_t , size_t, double)`| Allocates an `m`$\times$`n` matrix, with entries initialized to a value. Default value is `0.`.|
|`alg::Matrix(alg::Vector, size_t, size_t)`| Creates a matrix encapsulating a vector, with optional arguments. If the optional arguments are omitted, the matrix is shaped into column by default. |
|`alg::Matrix(alg::Vector, alg::Shape)`| Creates a matrix encapsulating a vector, but with shape defined by the `Shape` argument. |
|`alg::Matrix(std::vector<double>, size_t, size_t)`| Same as above. |
|`alg::Matrix(std::vector<double>, alg::Shape)`| Same as above. |
|`alg::Matrix(alg::Matrix)`| Copy constructor |
|`alg::Matrix(std::initializer_list<std::initializer_list<double>>)` | Allows initialization via initializer lists.|
|||

### Methods
| | |
|--|--|
|`at(size_t)`| Returns a reference to the element at given position in the flattened matrix (linewise).|
|`at(size_t, size_t)`| Returns a reference to the element at position $i,j$ for given arguments `i` and `j`.
|`begin()` | Returns an iterator pointing to the first element. The iterator runs linewise, that is, line by line. |
|`beginCol()`| Returns an iterator pointing to the first column.|
|`beginRow()`| Returns an iterator poiting to the first row.|
|`clear()` | Clears the contents of the object. Shape is set to `{ 0, 0 }`. |
|`col(size_t)`| Returns a Column object comprised of the specified column. This object efficiently gives references to entries in the matrix and can be used to modify the original matrix. |
|`concatenate(alg::Matrix)` | Concatenates the given matrix horizontally to this. |
|`det()` | Returns the determinant of this matrix if it's a square matrix. |
|`end()` | Returns an iterator pointing to just past the last element. The iterator runs linewise, that is, line by line. |
|`endCol()`| Returns an iterator pointing to just past the last column.|
|`endRow()`| Returns an iterator pointing to just past the last row.|
|`getInternalStdVector()`| Returns a reference to the underlying data container. |
|`getShape()`| Returns an `alg::Shape` struct with fields `m` and `n`, the number of lines and columns of the matrix.|
|`insert(alg::Iterator, double)`| Inserts a value at the position specified by the iterator. If the shape of the matrix was not planned to accommodate the insertion, it will have to be corrected with `setShape`.|
|`inv()` | Returns the inverse of this matrix.  |
|`isColumn()` | Returns true if the matrix has a single column, false otherwise. |
|`isScalar()` | Returns true if the matrix is $1\times1$, false otherwise. |
|`isRow()` | Returns true if the matrix has a single row, false otherwise. |
|`isVector()` | Returns true if the matrix has either a single column or a single row. |
|`ncols()`| The number of columns in the matrix. |
|`norm(double)`| The $p$-norm of the vector. Argument is the $p$ parameter. Default is `2.`.|
|`nrows()`| The number of rows in the matrix. |
|`push_back(alg::Vector)`| Appends a new line to the bottom of the matrix. |
|`reshape(alg::Shape)`| Returns a copy of the matrix with shape given by the argument.|
|`reshape(size_t, size_t)`| Returns a copy of the matrix with shape given by the two arguments.|
|`row(size_t)`| Returns a Row object containing the specified row. This object efficiently gives references to entries in the matrix and can be used to modify the original matrix. |
|`setShape(alg::Shape)`| Changes the shape of the matrix with shape given in the argument.|
|`setShape(size_t, size_t)`| Changes the shape of the matrix with the given parameters.|
|`size()`| Returns the total number of entries in the matrix. Same result as `nrows() * ncols()`.|
|`slice(size_t, size_t, size_t, size_t)` | Returns a submatrix of this matrix. Arguments are: upper row index, bottom row index, leftmost column index and rightmost column index. Bottommost and rightmost indexes are exclusive. |
|`t()`| Returns a transposed copy of the matrix.|
|`to_string()`| Returns a representation of the Matrix as `std::string`. |

### Functions on `alg::Matrix`
|||
|-|-|
|`alg::det(alg::Matrix)`| Returns the determinant of a matrix.|
|`alg::norm(alg::Matrix, double)` | Gives the $p$-norm of `m`, with the $p$ parameter given in the second argument. |
|`alg::t(alg::Matrix)` | Returns the transpose of a matrix. |

### Operators
Mathematical
+  `+` commutative. When used with:
	+ scalars
		+  Shifts every element by the scalar.  Return type: `alg::Matrix`
	+ vectors
		+ Sums only if matrix is row/column of same length as the vector. Return type: `alg::Vector`
	+ matrices
		+ Sums both matrices if they are of same shape.  Return type: `alg::Matrix`
+ `-` if
	+ unary:
		+ scales each element by -1.
	+ binary:
		+ same behavior as `+`.
+ `*` Only commutes with scalars. When used with:
	+ scalars
		+ Scales every element by the scalar. Return type: `alg::Matrix`
	+ vectors	
		+ If operands are of appropriate shapes, does product and returns the resulting vector (of same length as the one being operated upon). Return type: `alg::Vector`
	+ matrices
		+ Does matrix multiplication and returns a new matrix of appropriate shape.  Return type: `alg::Matrix`
+ `/` non-commutative. Can only be used with:
	+ scalars
		+ same behavior as `*`.
+ `+=` `-=` `*=` `/=` also avaliable. 

Access
+ `[]` Returns a Row object, the row at the respective index. 
+ `()` Performant way to get the behavior of `[][]`. Example: `matrix(i, j) == matrix[i][j]`.

Assignment
+ `=` When used with:
	+ matrices
		+ If shapes match, appropriately copies the contents of the right-hand side matrix to the left-hand side one.
	+ vectors
		+ If matrix is a either row or column and of same length, appropriately copies the contents of the right-hand side vector to the left-hand side matrix.
	+ row- or column-objects (special vectors)
		+ As special kinds of vector, the matrix is further restricted to be of the respective kind in order to proceed with the copy, not of just of either kind.

Miscellaneous
+ `<<` Appropriately sends matrix objects to the designated stream.

## Convenience functions
|||
|-|-|
|`alg::arange(double, double, double)`| Returns an `alg::Vector` generated by the arguments. First is the upper limit of the sequence, exclusive. Second is the lower limit of the sequence, inclusive. Last is the interval between elements (step size). Second and last arguments have defaults: `0.` and `1.`, respectively. Can be used to generate a matrix by being passed along with a shape parameter to the matrix constructor.| 
|`alg::I(size_t)`| Returns the identity matrix of given dimension.|

## Internal types: alg::Row and alg::Column

When the user needs to iterate over the rows or columns of the matrix, or simply needs access to a particular line of the matrix, special objects are created that hold references to all elements situated in that particular line, so that operations and modifications can be done in place with STL functions or user-defined routines. Implementation-wise the objects don't really keep data on each element of the line, they are just convenient devices to iterate over those elements in the matrix.

As such, these objects support iteration and the same numerical operations as `alg::Vector`. However, as to shape restrictions for those operations, they're as restricted as matrices. That is, while for instance an `alg::Vector` can be on either side of a dot product with an `alg::Matrix`, and `alg::Column` can only ever multiply a matrix from the left if the matrix is row-like. They also don't support insertions or removals. 

The most straightforward way to obtain such an object is by calling the methods `row` and `col` of `alg::Matrix`, such as in `auto row = matrix.row(3);`. Another "legal" manner is to iterate over the matrix with `beginRow`, `endRow`, `beginCol`, and `endCol`.

Following is a list of methods for these objects:
- `begin()`
- `end()`
- `getShape()` (returns the shape of the original matrix)
- `nrows()` (for `alg::Row`)
- `norm()` and `alg::norm(<row/column>)`
- `nrows()` (for `alg::Column`)
- `operator[]`
- `operator<<`
- `operator=` (with Vector, Matrix and same type)
- `size()`
- `to_string()`