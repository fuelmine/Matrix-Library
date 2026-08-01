#include <iostream>
#include "Matrix.hpp"

int main(){

    Matrix A = Matrix::identity(3);
    
    std::cout << A << "\n";
    
    double det = A.determinant();
    std::vector<double> sol = A.solve({1, 2, 3});
    std::cout << "The determiante of identity matrix is: " << det << "\n"; 
    for(int i = 0; i < sol.size(); ++i){
        std::cout << sol[i] << "\n";
    }

    return 0;
}