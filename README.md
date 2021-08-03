# Linear algebra in C++


The `linear_algebra.hh` header implements functionality and a numerical interface for linear algebra objects and constructs. Here is a short list of some of the things you can do with the functionality defined in this header:

- Easily declare matrix and vector objects.

````C++
#include "linear_algebra.hh"
#include<iostream>
using namespace::alg;
using std::cout;
int main() {
	Matrix<double> M = { {1, 2}, {3, 4} };
	cout << M;
}
````

O pacote `linear_algebra.hh` implementa 
funcionalidades de álgebra linear e de qualidade de vida para inerfaces numéricas. Como esse
é um tipo de projeto que pode ficar para a posteridade, ele foi todo implementado em inglês.
Segue uma lista das principais funcionalidades implementadas nesse pacote:
- Capacidade de *printar* vetores e matrizes
- Operações unárias de álgebra linear (transposição, norma, inversão de sinal etc)
- Operações binárias de álgebra linear com linguagem simples. A intenção é reproduzir
    a naturalidade de operações matriciais de pacotes como o Numpy ou o MATLAB. Isto é,
    o produto interno entre dois vetores **_v_** e **_u_** deve ser escrito 
    simplesmente utilizando o operador `*` utilizado para multiplicação: 
`double w = v * u;`
 - Utilidades. Exemplos: geração automática de objetos representando a matriz identidade
    *I<sub>n* para o argumento *n* ∈ ℤ<sup>+</sup>; rotina de
    inversão de matrizes etc;

> ###### Observação: 
> Todas as definições inclusas nesse pacote são paramétricas, de modo que o cabeçalho de definição do template (`template<typename T>`) será omitido neste documento.

Devido a possibilidade de o usuário encontrar utilidade em utilizar as definições desse
pacote e pela necessidade de documentação do mesmo, uma breve descrição será feita das 
funções e objetos mais importantes definidos.

> ###### Observação: 
> Todos os métodos e funções retornam uma nova cópia do objeto, deixando
os originais intactos.

- `Matrix<T>` Define um objeto representando uma matriz que se comporta quase exatamente
    como `vector<vector<T>>`, cujos elementos podem ser acessados pelo operador colchete 
    (`matriz[i][j]` ou `matriz[i]`) ou parêntese (`matriz(i, j)`). 
    
    Construtores:
    
     -   `{Matrix()`
     -   `{Matrix(std::vector<T> user_v, bool column = true)`
     -    `{Matrix(std::vector<std::vector<T>> user_m)`
     -    `{Matrix(Matrix<T>* user_m)`
     -    `Matrix(size_t m, size_t n)`
     
    Métodos notáveis da classe:
    
    -    `size_t n_rows()` A quantidade de linhas da matriz
    -    `Matrix<T> concatenate(Matrix<T> other)` Concatena duas matrizes
    -    `Matrix<T> slice(size_t, size_t, size_t, size_t)]` Retorna uma submatriz
    -    `Matrix<T> inverse()` Retorna a matriz inversa
    -    `Matrix<T> t(Matrix<T> M)` Retorna a matriz transposta
    -    `double norm()` Retorna a norma caso a matriz seja linha ou coluna 
    -    `static double norm(std::vector<T> v)` Sobrecarga do método anterior
    
    
   - `Matrix<T> I(size_t m)` Retorna uma matriz identidade de dimensão _m_
   - `std::vector<T> range(unsigned int n)` Retorna um vetor _{1,...,n}_
   - `Matrix<T> t(std::vector<T> v)` Retorna uma matriz coluna a partir do vetor **_v_**
   -  **Sobrecarga de operadores** Operadores sobrecarregados: `+ - * += -= *= <<`

As outras definições importantes dizem respeito à sobrecarga dos operadores aritméticos 
`+ - *`. No caso dos operadores adição, subtração e sinal negativo (existe uma 
diferença entre os dois últimos em programação), todas as operações se comportam como o
esperado para álgebra linear. O operador `<<` emite 
corretamente objetos do tipo`vector<T>` e \`Matrix<T>` para a _stream_
designada.
