# Lincpp

### Linear algebra in C++

## Summary
- [Lincpp](#lincpp)
		- [Linear algebra in C++](#linear-algebra-in-c)
	- [Summary](#summary)
	- [Introduction](#introduction)
			- [Example of use](#example-of-use)
	- [Description](#description)
	- [The `Vector` class](#the-vector-class)
		- [Constructors](#constructors)
		- [Methods](#methods)
		- [Functions on `alg::Vector`](#functions-on-algvector)
		- [Operators](#operators)
	- [The `Matrix` class](#the-matrix-class)
		- [Shape](#shape)
		- [Constructors](#constructors-1)
		- [Methods](#methods-1)
		- [Functions on `alg::Matrix`](#functions-on-algmatrix)
		- [Operators](#operators-1)
	- [Convenience functions](#convenience-functions)
	- [Internal types: alg::Row and alg::Column](#internal-types-algrow-and-algcolumn)


## Introduction
This library implements functionality and a numerical interface for linear algebra objects and constructs, with aim to make numerical computations in C++ as straightforward to work with as in popular environments such as MATLAB and Numpy. 

The main goal is ease of use, readability and numerical efficiency. Since these are somewhat contradictory goals, a compromise has been made. In order to preserve a straightforward mathematical syntax, the "return by value/copy" approach has been used. This way, a complex mathematical expression can be written on the right hand side of an attribution by exploiting operator overloads for algebraic classes.

The idea is to allow for rapid prototyping of numerical computations à la Python and, if and when necessary, the developer can easily locate and rewrite the mission-critical bottlenecks by hand without giving up on the very useful numeric data structures.

## Examples

Below is a short list of examples of what you can do with the functionality defined here:

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
std::cout << A[0].norm() << "\n"; // Note that this is the norm of row 0 of A
std::cout << A.inverse();         // [TODO:] You can also take `A.norm(2)`
````
Output:
````
2.23607

{{ -2 1 }
 { 1.5 -0.5 }}

````

- Operate on rows and columns by reference.

````C++
#include <numeric>

Matrix B(3, 3);
std::iota(B.begin(), B.end(), 1);

std::for_each(B.beginCol(), B.endCol(), [](auto column) {
	column += Vector(3, 10);
	column[1] += 100;
});

std::cout << B;
````
Output:
````
{{ 11 12 13 }
 { 114 115 116 }
 { 17 18 19 }}
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
### A more complex example
For a full example of application, check [this submodule](https://github.com/CoralBleaching/Non-linear-optimization) on non-linear optimization. It makes use of these definitions on its core routines.

Below is one function defined in that project. It takes a functional object `G` that represents a vector function, for instance, the gradient of a multivariable function, $G(x) = \nabla f(x)$, and returns another function that computes its Jacobian at any specific point via a central difference formula:

$$
H_{ij} = \frac{[G(\vec x + eI_j) - G(\vec x - eI_j)]_i + [G(\vec x + eI_i) - G(\vec x - eI_i)]_j}{4e},
$$

where $H = (H_{ij})$ is the resulting Jacobian matrix, $I$ is the identity matrix, and $e$ is the desired precision. Here we call the Jacobian $H$ instead of $J$ because in this project the function is tasked to find the Hessian of a function $f(x)$ from its gradient (or its second derivative from its first derivative, in case of a $x$ being a scalar).

````C++
#include <functional>
// notice how the expression within the innermost loop automatically evaluates to a double
std::function<Matrix(Vector)> jacobian(std::function<Vector(Vector)> G, int n, double e = 1e-4)
{
    return [G, n, e](Vector x) {
        Matrix H(n, n);
        Matrix M = e * I(n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                H(i, j) = (G(x + M[j]) - G(x - M[j]))[i] / (4 * e) + (G(x + M[i]) - G(x - M[i]))[j] / (4 * e);
        H = 0.5 * (H + t(H)); // this could be more efficiently but less conveniently done by hand
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
> Nevertheless, there are instances when you don't want the result of a dot product to be reduced to an object of type `alg::Vector`, or the result might still be a matrix, for example: $(v_1, \dots, v_m)^T\cdot(u_1,\dots,u_m)$. For this reason, operations between two objects of type `alg::Matrix` still always return a new object of type `alg::Matrix`. Operations of the type `A.col(j) * A.row(i)` (where `A` is necessarily a square `alg::Matrix`) also yield an `alg::Matrix` return type.

**The header defines its own namespace, `alg`.**

Vector properties and behavior were encapsulated in a class called `Vector`. Matrix properties and behavior were encapsulated in a class called `Matrix`.
> ##### Note
> Methods and operators return new copies of the objects. Few methods act in-place, mostly insertion and access methods. For instance, the `reshape` method returns a copy, but the `setShape` method acts in-place.

Vectors do not distinguish between column or row types. When coupled with other objects like matrices and other vectors, the appropriate behavior is inferred. For example, both `A*v` and `v*A`, where `A` is a `Matrix` and `v` is a `Vector`, give the expected results, as `v` is conceptually transposed accordingly.

## The `Vector` class

`alg::Vector` defines an object representing a mathematical vector which behaves in C++ similarly to an 1-D NDarray in Python. Its elements can be accessed by (const or non-const) reference via brackets or `at` method (`vector[i]` or `vector.at(i)`).

### Constructors
|||
|-|-|
| `alg::Vector(size_t n, double f)` | Default constructor (defaults: `0`, `0.`). Explicit. Allocates a vector of length `n` and initializes all entries with the same value `f`. |
| `alg::Vector(std::vector<double>)` | Conversion constructor for `std::vector`. Explicit. |
| `alg::Vector(alg::Matrix)` | Flattens the Matrix, losing shape information. Explicit. |
| `alg::Vector(Iterator, Iterator)` | Initializes a new `alg::Vector` with values from the range defined by the iterators. Iterators must be the ones built-in for `Vector`, `Matrix`, rows or columns. |
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
|`push_back(double)`| Appends a new element after the last element. Note that this changes the mathematical properties of the object. |
|`row()`| Returns a row `alg::Matrix` with the same values. Same as `t()`. |
|`size()`| Returns the length of the vector. |
|`slice(size_t s, size_t f)` | Returns a vector comprised of a subset of elements within the range `s`:`f` of the original vector. Will throw if limit exceeds the vector's lenght. |
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
+ `<<` Sends a representation of the vector object to the designated stream.

## The `Matrix` class

`alg::Matrix` defines an object representing a matrix in row-major order (which behaves like an array of vector objects although it's implemented with contiguous data). Its elements can be accessed via parentheses (`matrix(i, j)`).  Brackets (`matrix[i][j]`) notation is possible techincally, but the parentheses method is the direct one. `matrix[i]` will return a `Row` object that can be further indexed to extract an element at position `j`.

### Shape

`alg::Shape` is a struct keeping a pair of fields `m`, `n` of the number of rows and columns, respectively, of a matrix. It's shorthand used for passing this information packaged together. Like any simple object it can be instantiated with brackets (initializer lists), e.g. `{2, 3}`.

It also features a comparison operator `==` for testing equality of shapes.

### Constructors
| | |
|--|--|
|`alg::Matrix(alg::Shape s, double f)`| Default constructor. Builds a matrix with shape `s` where each element is `f`. Defaults are `{ 0, 0 }` and `0.`, respectively. |
|`alg::Matrix(size_t m, size_t n, double f)`| Allocates an `m`$\times$`n` matrix, with entries initialized to a value `f`. Default value is `0.`.|
|`alg::Matrix(alg::Vector v, size_t m, size_t n)`| Creates a matrix encapsulating a vector, with optional arguments. If the optional arguments are omitted, the matrix is shaped into column by default. |
|`alg::Matrix(alg::Vector, alg::Shape)`| Creates a matrix encapsulating a vector, but with shape defined by the `Shape` argument. |
|`alg::Matrix(std::vector<double>, size_t, size_t)`| Same as above. |
|`alg::Matrix(std::vector<double>, alg::Shape)`| Same as above. |
|`alg::Matrix(alg::Matrix)`| Copy constructor |
|`alg::Matrix(std::initializer_list< std::initializer_list<double> >)` | Allows construction from initializer lists.|
|||

### Methods
| | |
|--|--|
|`at(size_t)`| Returns a reference to the element at given position in the flattened matrix (row-major order).|
|`at(size_t i, size_t j)`| Returns a reference to the element at row `i`, column `j`.
|`begin()` | Returns an iterator pointing to the first element. The iterator runs in row-major order, that is, line by line. |
|`beginCol()`| Returns an iterator pointing to the first column (the entire column).|
|`beginRow()`| Returns an iterator poiting to the first row (the entire row).|
|`clear()` | Clears the contents of the object. Shape is set to `{ 0, 0 }`. |
|`col(size_t)`| Returns a Column object comprised of the specified column. This object efficiently gives references to entries in the matrix and can be used to modify the original matrix. |
|`concatenate(alg::Matrix)` | Generates a new Matrix resulting from concatenating the given matrix horizontally to this. |
|`det()` | Returns the determinant of this matrix if it's a square matrix. |
|`end()` | Returns an iterator pointing to just past the last element. The iterator runs in row-major order, that is, line by line. |
|`endCol()`| Returns an iterator pointing to just past the last column.|
|`endRow()`| Returns an iterator pointing to just past the last row.|
|`getInternalStdVector()`| Returns a const reference to the underlying data container. |
|`getShape()`| Returns an `alg::Shape` struct with fields `m` and `n`, the number of lines and columns of the matrix.|
|`length()`| Returns the total number of entries in the matrix. Same result as `nrows() * ncols()`.|
|`insertColumn(iterator b, iterator e, size_t j = <end>)`| Inserts an entire column defined by the range [`b`, `e`) at the position specified by `j`. Default position is after the last column. If the matrix is empty, it will insert and reshape the matrix accordingly. Will throw if the length is incompatible. |
|`insertColumn(<Vector, Matrix, Row or Column> c, size_t j = <end>)`| Inserts an entire array `c` as a column at the position specified by `j`. Default position is after the last column. If the matrix is empty, it will insert and reshape the matrix accordingly. Will throw if the length is incompatible. |
|`insertRow(iterator b, iterator e, size_t i = <end>)`| Inserts an entire row defined by the range [`b`, `e`) at the position specified by `i`. Default position is after the last row. If the matrix is empty, it will insert and reshape the matrix accordingly. Will throw if the length is incompatible. |
|`insertRow(<Vector, Matrix, Row or Column> r, size_t i = <end>)`| Inserts an entire array `r` as a row at the position specified by `i`. Default position is after the last row. If the matrix is empty, it will insert and reshape the matrix accordingly. Will throw if the length is incompatible. |
|`inv()` | Returns the inverse of this matrix.  |
|`isColumn()` | Returns true if the matrix has a single column, false otherwise. |
|`isRow()` | Returns true if the matrix has a single row, false otherwise. |
|`isScalar()` | Returns true if the matrix is $1\times1$, false otherwise. |
|`isVector()` | Returns true if the matrix has either a single column or a single row. |
|`ncols()`| The number of columns in the matrix. |
|`norm(double)`| The $p$-norm of the vector. Argument is the $p$ parameter. Default is `2.`.|
|`nrows()`| The number of rows in the matrix. |
|`reshape(alg::Shape)`| Returns a copy of the matrix rearranged with shape given by the argument.|
|`reshape(size_t m, size_t n)`| Returns a copy of the matrix rearranged to `m` rows and `n` columns. |
|`row(size_t)`| Returns a Row object containing the specified row. This object efficiently gives references to entries in the matrix and can be used to modify the original matrix. |
|`setShape(alg::Shape)`| Changes the shape of the matrix with shape given in the argument.|
|`setShape(size_t m, size_t n)`| Changes the shape of the matrix to `m` rows and `n` columns. |
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