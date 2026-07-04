#pragma 0
#include <vector>
#include <initializer_list>

class Matrix{
private:
    size_t rows_;
    size_t cols_;
    std::vector<double> matrix;
    //Helper Function that gives us the index of (i,j) Note: 0-indexed
    size_t index(size_t row, size_t column) const;
public:
    //Constructor & Destructor
    Matrix(size_t row, size_t column) : rows_{row}, cols_{column}, matrix(row * column, 0){ }
    Matrix(size_t row, size_t column, double value) : rows_{row}, cols_{column}, matrix(row * column, value) { }
    Matrix(size_t row, size_t column, std::initializer_list<double> values) : rows_{row}, cols_{column}, matrix(values){
        if(values.size() != row * column){
            throw std::invalid_argument("Wrong number of values");
        }
    }

    //operator overloading
    double operator()(size_t i, size_t j) const;
    double& operator()(size_t i, size_t j);

    //include matrix arithmetic
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const; //O(n^3) time complexity 
    Matrix operator*(double scalar) const;

    //Linear Algerbra Functions
    bool isSquare() const;
    Matrix transpose() const;
    double trace() const;
    double determinate() const;
    size_t rank() const;
    std::vector<double> solve(const std::vector<double> &b) const; //O(n^3) time complexity 

    //print function for matrices
    void print() const;

    //getters for the size of row and columns
    size_t rows() const{return rows_;}
    size_t cols() const{return cols_;}

};

Matrix operator*(double scalar, const Matrix &other);