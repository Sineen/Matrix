//
// Created by sin_een on 1/14/19.
//

#ifndef EX3_MATRIX_HPP
#define EX3_MATRIX_HPP

#include <iostream>
#include <vector>

using namespace std;

template <class T>
class Matrix
{
private:
    unsigned int columns;
    unsigned int rowNumbers;
    unsigned int size;
    vector<T> matrix;

public:

    /**
     * empty constructor that building a matrix of one
     */
    Matrix(): columns(1), rowNumbers(1), size(1), matrix(0){}

    Matrix(unsigned int rows, unsigned int cols);

    Matrix(vector<T> other);

    Matrix(unsigned int rows, unsigned int cols, const vector<T>& cells);

    ~Matrix();

    Matrix<T>& operator=(const Matrix<T> &other);

    /**
     * subtracts two matrixes  and puts the vlaues  on teh matrix in the class
     * @param other a second matrix
     * @return  matrix that is sum of both matrixes
     */
    Matrix<T> operator+(const Matrix<T> &other);

    /**
     * multiplications two matrixes  and puts the vlaues  on teh matrix in the class
     * @param other a second matrix
     * @return matrix that is the subtraction between th e2 matrixes
     */
    Matrix<T> operator-(const Matrix<T> &other);

    /**
     * addes the object of this class to other adn makes a new matrix and returns that
     * @param other another matrix3D
     */
    Matrix<T> operator*(const Matrix &other);

    /**
     *
     */
    Matrix<T> operator==(const Matrix<T>& other);

    /**
     *
     */
    Matrix<T> operator!=(const Matrix<T>& other);

    /**
     * only works for square matrixes ( colum = row) other wise retrun exception
     * @return a new matrix that isi trasposed to our matrix
     */
    Matrix<T> trans();


    /**
     * @return true if the matrix is square
     */
    bool isSquareMatrix();


    /**
     *
     */
     friend Matrix<T> operator<<(ostream& oc, Matrix<T> matrix1);

     /**
      *
      */
     <T> operator (unsigned int, unsigned int);

     /**
      *
      */
     const <T> operator()(unsigned int, unsigned int);


/// \returns the begin iterator
    typename std::vector<T>::iterator begin() {
        return matrix.begin();
    }

/// \returns the end iterator
    typename std::vector<T>::iterator end() {
        return matrix.end();
    }

    unsigned int rows();

    unsigned int cols();

};

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols)
{
    this->size = rows*cols;
    this->rowNumbers = rows;
    this->columns = cols;
    this->matrix.clear();
    for(int i = 0; i < size; i++ )
    {
        matrix.push_back(T());
    }
}

template<class T>
Matrix<T>::Matrix(vector<T> other)
{
    this->size = other.size;
    this->rowNumbers = other.rowNumbers;
    this->columns = other.columns;
    this->matrix.clear();
    for(int i = 0; i < size; i++ )
    {
        matrix.push_back(T());
    }
}


#endif //EX3_MATRIX_HPP
