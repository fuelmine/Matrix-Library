/*
    Project Name: Matrix Library
    Descripton: Implementation of a matrix library
*/
#include <iostream>
#include <vector>
#include <exception>
#include <cmath>
#include <algorithm>
#include "Matrix.hpp"


size_t Matrix::index(size_t row, size_t column) const{
    if(rows_ <= row || cols_ <= column)
        throw std::out_of_range("Matrix index out of range");
    
    return (row * cols_) + column;
}

Matrix Matrix::identity(size_t n){
    Matrix I(n, n);
    for(size_t i{}; i < n; ++i){
        I(i, i) = 1.0;
    }
    return I;
}
//All the operator overloading.
//1. This return the element at (i,j)
const double& Matrix::operator()(size_t i, size_t j) const{
    return matrix_[index(i, j)];
}
//2. This allow us to put a value at element (i, j)
double& Matrix::operator()(size_t i, size_t j){
    return matrix_[index(i, j)];
}
//3. The sum of matrices A + B
Matrix Matrix::operator+(const Matrix &other) const{
    //check that the two matrices have the same size
    if(rows_ != other.rows() || cols_ != other.cols())
        throw std::invalid_argument("Matrix dimensions must match");
    //now sum each the matrices
    Matrix sum(rows_, cols_);
    for(size_t i{}; i < rows_; ++i){
        for(size_t j{}; j < cols_; ++j){
            sum(i, j) = (*this)(i,j) + other(i, j);
        }
    }
    return sum;
}
//4. The difference of matrices A - B
Matrix Matrix::operator-(const Matrix &other) const{
        //check that the two matrices have the same size
    if(rows_ != other.rows() || cols_ != other.cols())
        throw std::invalid_argument("Matrix dimensions must match");
    //now sum each the matrices
    Matrix diff(rows_, cols_);

    for(size_t i{}; i < rows_; ++i){
        for(size_t j{}; j < cols_; ++j){
            diff(i, j) = matrix_[index(i,j)] - other(i, j);
        }
    }
    return diff;
}
//5. The product of matrices AB
Matrix Matrix::operator*(const Matrix &other) const{
    if(cols_ != other.rows())
        throw std::invalid_argument("Matrix dimensions must match");
    
    Matrix result(rows_, other.cols(), 0.0);
    //get product O(n^3)
    for (size_t i{}; i < rows_; ++i) {
        for (size_t k{}; k < cols_; ++k) {
            for (size_t j{}; j < other.cols_; ++j) {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}
//6. The product of a matrix and scalar on the right side Ma
Matrix Matrix::operator*(double scalar) const{
    Matrix result(rows_, cols_);
    for(size_t i{}; i < rows_; ++i){
        for(size_t j{}; j < cols_; ++j){
            result(i, j) = scalar * (*this)(i , j);
        }
    }
    return result;
}
//7. The product of a matrix and scalr on the left side aM
Matrix operator*(double scalar, const Matrix &other){
    return other * scalar;
}

void Matrix::print() const{
    for(size_t i{}; i < rows_; ++i){
        for(size_t j{}; j < cols_; ++j){
            std::cout << (*this)(i,j) << " ";
        }
        std::cout << "\n";
    }
}

//returns true if the a matrix is square
bool Matrix::isSquare() const{
    if(rows_ != cols_) return false;

    return true;
}

//This function tranposes the Matrix
Matrix Matrix::transpose() const{
    Matrix transposed(cols_, rows_);
    for(size_t i{}; i < cols_; ++i){
        for(size_t j{}; j < rows_; ++j){
            transposed(i, j) = (*this)(j, i);
        }
    }
    return transposed;
}

//This function returns the trace of a matrix
double Matrix::trace() const{
    if(!isSquare()){
        throw std::invalid_argument("Matrix must be square");
    }
    double result{0.0};
    for(size_t i{}; i < rows_; ++i){
        result += (*this)(i, i);
    }
    return result;
}

//This function gives the determinate of a matrix using the upper triangle matrix.
double Matrix::determinant() const{
    if(!isSquare()){
        throw std::invalid_argument("Matrix must be square");
    }
    Matrix A = *this;
    size_t n = rows_;
    size_t numSwaps{};

    for(size_t i{}; i < n; ++i){
        //find the appropraite pivot
        size_t pivot = i;
        for(size_t j{i + 1}; j < n; ++j){
            if(std::abs(A(j, i)) > std::abs(A(pivot, i))) pivot = j;
        }
        if(std::abs(A(pivot, i)) < 1e-12){
            return 0.0;
        }
        //swap rows
        if(pivot != i){
            for(size_t k{i}; k < n; ++k){
                std::swap(A(i, k), A(pivot, k));
            }
            numSwaps++;
        }
        //eliminate
        for(size_t j{i + 1}; j < n; ++j){
            double factor{A(j, i)/ A(i, i)};
            for(size_t k{i}; k < n; ++k){
                A(j, k) -= factor * A(i, k);
            }
        }
    }
    //Get the product of diagonal
    double det{1.0};
    for(size_t i{}; i < n; ++i){
        det *= A(i,i);
    }

    //if numSwap is odd return the negative
    if(numSwaps % 2) return -det;

    return det;
}

// Gaussian Elimiation Function. Returns solution but if matrix is singular throws exception.
std::vector<double> Matrix::solve(const std::vector<double> &b) const{
    if(!isSquare()){
        throw std::invalid_argument("Matrix must be square");
    }

    if(b.size() != rows_){
        throw std::invalid_argument("Invalid b size");
    }

    Matrix A = *this;
    std::vector<double> rhs = b;
    size_t n = rows_;

    for(size_t i = 0; i < n; ++i){
        //find pivot 
        size_t pivot = i;
        for(size_t j = i + 1; j < n; ++j){
            if(std::abs(A(j, i)) > std::abs(A(pivot, i))) pivot = j;
        }
        //check if A is singular matrix 
        if(std::abs(A(pivot, i)) < 1e-12){
            throw std::runtime_error("Matrix is singular");
        }
        //swap rows
        if(pivot != i){
            for(size_t k = i; k < n; ++k){
                std::swap(A(i, k), A(pivot, k));
            }
            std::swap(rhs[i], rhs[pivot]);
        }
        //Eliminate the variable 
        for(size_t j = i + 1; j < n; ++j){
            double factor = A(j, i) / A(i, i);
            for(size_t k{i}; k < n; ++k){
                A(j, k) -= factor * A(i, k);
            }
            rhs[j] -= factor * rhs[i];
        }
    }

    std::vector<double> x(n, 0.0);
    //Substitution 
    for(size_t i{n}; i-- > 0; ){
        double temp{0.0};
        for(size_t j{i + 1}; j < n  ; ++j) temp += A(i, j) * x[j];

        x[i] = (rhs[i] - temp) / A(i, i);
    }

    return x;
}

size_t Matrix::rank() const{
    Matrix A = *this;
    size_t rank{};

    for(size_t col{}; col < cols_ && rank < rows_; ++col){
        size_t pivot = rank;

        for(size_t row{rank + 1}; row < rows_; ++row){
            if(std::abs(A(row, col)) > std::abs(A(pivot, col))) pivot = row;
        }

        if(std::abs(A(pivot, col)) < 1e-12) continue;

        if(pivot != rank){
            for(size_t j{col}; j < cols_; ++j){
                std::swap(A(rank, j), A(pivot, j));
            }
        }
        for(size_t row = rank + 1; row < rows_; ++row){
            double factor = A(row, col) / A(rank, col);
            for (size_t j = col; j < cols_; ++j) {
                A(row, j) -= factor * A(rank, j);
            }
        }

        ++rank;
    }
    return rank;
}


int main(){

    Matrix A = Matrix::identity(3);
    A.print();
    
    double det = A.determinant();
    
    std::vector<double> sol = A.solve({1, 2, 3});
    std::cout << "The determiante of identity matrix is: " << det << "\n"; 
    for(int i = 0; i < sol.size(); ++i){
        std::cout << sol[i] << "\n";
    }

    return 0;
}
