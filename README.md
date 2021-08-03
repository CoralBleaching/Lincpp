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
|`size_t n_rows()`| A quantidade de linhas da matriz|
|`Matrix<T> concatenate(Matrix<T>)` | Concatena duas matrizes |
|`Matrix<T> slice(size_t, size_t, size_t, size_t)]` | Retorna uma submatriz |
|`Matrix<T> inverse()` | Retorna a matriz inversa |
|`Matrix<T> t(Matrix<T> M)` | Retorna a matriz transposta |
|`double norm()` | Retorna a norma caso a matriz seja linha ou coluna |
|`static double norm(std::vector<T> v)` |Sobrecarga do método anterior|
|`Matrix<T> I(size_t m)`| Retorna uma matriz identidade de dimensão _m_|
| `std::vector<T> range(unsigned int n)` | Retorna um vetor _{1,...,n}_|
|`Matrix<T> t(std::vector<T> v)`| Retorna uma matriz coluna a partir do vetor **_v_**|

#### Sobrecarga de operadores

Operadores sobrecarregados: `+ - * += -= *= <<`

As outras definições importantes dizem respeito à sobrecarga dos operadores aritméticos 
`+ - *`. No caso dos operadores adição, subtração e sinal negativo (existe uma 
diferença entre os dois últimos em programação), todas as operações se comportam como o
esperado para álgebra linear. O operador `<<` emite 
corretamente objetos do tipo`vector<T>` e \`Matrix<T>` para a _stream_
designada.
